# Created by alexandra at 14/08/2025
"""
Paired comparison of two classification models evaluated on the same samples.

Inputs (CSV):
  - regression_results.csv  (columns: PatientId, ichorDetection, Predicted_label OR Predicted_class)
  - alfassay_results.csv    (columns: PatientId, ichorDetection, Predicted_label OR Predicted_class)

What it does:
  1) Inner-merge by PatientId.
  2) Uses ichorDetection as the binary ground truth (0/1).
  3) Runs paired tests:
       - Accuracy, Sensitivity, Specificity: McNemar exact + paired bootstrap CI/p-value
       - PPV, NPV: paired bootstrap CI/p-value (McNemar not applicable)
  4) Writes results to paired_comparison_stats.csv with BH FDR-adjusted p-values.

Notes:
  - If both files contain ichorDetection and they disagree for some PatientId, the script warns
    and uses the ALFAssay value as the truth (change easily if you prefer the regression one).
  - For AUC comparison, you’d need continuous scores; use DeLong’s test (not included here).
"""

import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import binomtest

# --------------------
# Config
# --------------------
REG_PATH  = "/Users/alexandra/PhD/PyCharmProjects/ALFAssayNN/Rscripts/plots_t35/validation_regression_predictions_with_classes.csv"   # <-- change to your file path
ALFA_PATH = "/Users/alexandra/PhD/PyCharmProjects/ALFAssayNN/Rscripts/plots_t35/validation_predictions_with_classes.csv"     # <-- change to your file path

OUT_CSV   = "/Users/alexandra/PhD/PyCharmProjects/ALFAssayNN/Rscripts/plots_t35/paired_comparison_stats.csv"
N_BOOT    = 10000
SEED      = 123

# --------------------
# Helpers
# --------------------
def pick_pred_col(df):
    """Return the name of the predicted CLASS column, trying common options."""
    candidates = ["Predicted_label", "Predicted_class", "pred_class", "y_pred"]
    for c in candidates:
        if c in df.columns:
            return c
    raise KeyError(f"None of {candidates} found in columns: {list(df.columns)}")

def metrics(y_true, y_pred):
    """Basic classification metrics from confusion counts."""
    y_true = y_true.astype(int)
    y_pred = y_pred.astype(int)
    tp = int(((y_true == 1) & (y_pred == 1)).sum())
    tn = int(((y_true == 0) & (y_pred == 0)).sum())
    fp = int(((y_true == 0) & (y_pred == 1)).sum())
    fn = int(((y_true == 1) & (y_pred == 0)).sum())
    acc  = (tp + tn) / (tp + tn + fp + fn) if (tp + tn + fp + fn) > 0 else np.nan
    sens = tp / (tp + fn) if (tp + fn) > 0 else np.nan
    spec = tn / (tn + fp) if (tn + fp) > 0 else np.nan
    ppv  = tp / (tp + fp) if (tp + fp) > 0 else np.nan
    npv  = tn / (tn + fn) if (tn + fn) > 0 else np.nan
    return {"TP": tp, "TN": tn, "FP": fp, "FN": fn,
            "Accuracy": acc, "Sensitivity": sens, "Specificity": spec, "PPV": ppv, "NPV": npv}

def mcnemar_exact_pairs(x_hit, y_hit):
    """
    x_hit, y_hit = binary correctness indicators in the *same* set of items.
    Returns two-sided exact McNemar p-value and discordant counts b,c.
    """
    b = int(((x_hit == 1) & (y_hit == 0)).sum())  # x correct, y wrong
    c = int(((x_hit == 0) & (y_hit == 1)).sum())  # x wrong, y correct
    n_disc = b + c
    if n_disc == 0:
        return 1.0, b, c
    p = binomtest(min(b, c), n_disc, 0.5, alternative="two-sided").pvalue
    return float(p), b, c

def paired_boot_ci(y_true, y1, y2, metric_fn, B=10000, seed=123):
    """
    Paired bootstrap CI and p-value for difference metric(y2) - metric(y1).
    metric_fn must accept (y_true, y_pred) and return a scalar metric.
    """
    rng = np.random.default_rng(seed)
    n = len(y_true)
    diffs = np.empty(B, dtype=float)
    for b in range(B):
        samp = rng.integers(0, n, size=n)
        yt  = y_true[samp]
        y1b = y1[samp]
        y2b = y2[samp]
        m1 = metric_fn(yt, y1b)
        m2 = metric_fn(yt, y2b)
        diffs[b] = m2 - m1
    diffs = diffs[~np.isnan(diffs)]
    lo, hi = np.percentile(diffs, [2.5, 97.5])
    p_boot = 2 * min((diffs <= 0).mean(), (diffs >= 0).mean())  # two-sided
    return float(lo), float(hi), float(p_boot)

# Metric functions for bootstrap
def acc_fn(y, yhat):  return (y == yhat).mean()
def sens_fn(y, yhat):
    mask = (y == 1)
    return (yhat[mask] == 1).mean() if mask.any() else np.nan
def spec_fn(y, yhat):
    mask = (y == 0)
    return (yhat[mask] == 0).mean() if mask.any() else np.nan
def ppv_fn(y, yhat):
    mask = (yhat == 1)
    return (y[mask] == 1).mean() if mask.any() else np.nan
def npv_fn(y, yhat):
    mask = (yhat == 0)
    return (y[mask] == 0).mean() if mask.any() else np.nan

def bh_fdr(pvals):
    """
    Benjamini–Hochberg FDR adjustment (returns array of adjusted p-values).
    NaNs are preserved in position.
    """
    p = np.array(pvals, dtype=float)
    n = np.sum(~np.isnan(p))
    adj = np.full_like(p, np.nan, dtype=float)
    if n == 0:
        return adj
    # ranks among non-NaN
    idx = np.where(~np.isnan(p))[0]
    p_non = p[idx]
    order = np.argsort(p_non)
    ranks = np.empty_like(order)
    ranks[order] = np.arange(1, n+1)
    # BH: p_i * m / rank
    adj_non = p_non * n / ranks
    # monotone
    adj_non = np.minimum.accumulate(adj_non[::-1])[::-1]
    adj_non = np.clip(adj_non, 0.0, 1.0)
    adj[idx] = adj_non
    return adj

# --------------------
# Load & merge
# --------------------
reg = pd.read_csv(REG_PATH)
alfa = pd.read_csv(ALFA_PATH)

pred_reg_col  = pick_pred_col(reg)
pred_alfa_col = pick_pred_col(alfa)

need = ["PatientId", "ichorDetection", pred_reg_col]
reg = reg.loc[:, need].rename(columns={"ichorDetection": "true_reg", pred_reg_col: "pred_reg"})

need = ["PatientId", "ichorDetection", pred_alfa_col]
alfa = alfa.loc[:, need].rename(columns={"ichorDetection": "true_alfa", pred_alfa_col: "pred_alfa"})

df = reg.merge(alfa, on="PatientId", how="inner")

# Resolve truth if two files disagree
if not np.array_equal(df["true_reg"].values, df["true_alfa"].values):
    mismatches = int((df["true_reg"] != df["true_alfa"]).sum())
    print(f"[WARN] {mismatches} samples have differing 'ichorDetection' between files; using ALFA values.")
    df["true"] = df["true_alfa"]
else:
    df["true"] = df["true_reg"]

df = df[["PatientId", "true", "pred_reg", "pred_alfa"]].copy()

# Coerce to ints and drop NAs
for c in ["true", "pred_reg", "pred_alfa"]:
    df[c] = pd.to_numeric(df[c], errors="coerce")
df = df.dropna(subset=["true", "pred_reg", "pred_alfa"]).astype({"true": int, "pred_reg": int, "pred_alfa": int})

n   = len(df)
pos = int((df["true"] == 1).sum())
neg = int((df["true"] == 0).sum())
print(f"n={n}, positives={pos}, negatives={neg}")

y = df["true"].to_numpy()
r = df["pred_reg"].to_numpy()
a = df["pred_alfa"].to_numpy()

# --------------------
# Per-model metrics
# --------------------
m_reg  = metrics(y, r)
m_alfa = metrics(y, a)

# --------------------
# Paired exact tests (McNemar)
# --------------------
# Accuracy (correct vs incorrect)
acc_reg_hit  = (r == y).astype(int)
acc_alfa_hit = (a == y).astype(int)
p_acc, b_acc, c_acc = mcnemar_exact_pairs(acc_reg_hit, acc_alfa_hit)

# Sensitivity (in positives only, hit = predict 1)
pos_mask = (y == 1)
sens_reg_hit  = (r[pos_mask] == 1).astype(int)
sens_alfa_hit = (a[pos_mask] == 1).astype(int)
p_sens, b_sens, c_sens = mcnemar_exact_pairs(sens_reg_hit, sens_alfa_hit)

# Specificity (in negatives only, hit = predict 0)
neg_mask = (y == 0)
spec_reg_hit  = (r[neg_mask] == 0).astype(int)
spec_alfa_hit = (a[neg_mask] == 0).astype(int)
p_spec, b_spec, c_spec = mcnemar_exact_pairs(spec_reg_hit, spec_alfa_hit)

# --------------------
# Paired bootstrap CIs (differences ALFA - REG)
# --------------------
acc_ci_lo,  acc_ci_hi,  acc_p_boot  = paired_boot_ci(y, r, a, acc_fn,  B=N_BOOT, seed=SEED)
sens_ci_lo, sens_ci_hi, sens_p_boot = paired_boot_ci(y, r, a, sens_fn, B=N_BOOT, seed=SEED)
spec_ci_lo, spec_ci_hi, spec_p_boot = paired_boot_ci(y, r, a, spec_fn, B=N_BOOT, seed=SEED)
ppv_ci_lo,  ppv_ci_hi,  ppv_p_boot  = paired_boot_ci(y, r, a, ppv_fn,  B=N_BOOT, seed=SEED)
npv_ci_lo,  npv_ci_hi,  npv_p_boot  = paired_boot_ci(y, r, a, npv_fn,  B=N_BOOT, seed=SEED)

# --------------------
# Assemble results
# --------------------
def make_row(name, reg_val, alfa_val, ci_lo, ci_hi, p_exact, p_boot):
    diff = alfa_val - reg_val
    return {
        "Metric": name,
        "Regression": round(reg_val, 4),
        "ALFAssay": round(alfa_val, 4),
        "Difference (ALFA-REG)": round(diff, 4),
        "95% CI of Difference": f"[{ci_lo:.4f}, {ci_hi:.4f}]",
        "Exact p (McNemar)": None if p_exact is None else float(p_exact),
        "Bootstrap p": float(p_boot)
    }

rows = [
    make_row("Accuracy",    m_reg["Accuracy"],    m_alfa["Accuracy"],    acc_ci_lo,  acc_ci_hi,  p_acc,  acc_p_boot),
    make_row("Sensitivity", m_reg["Sensitivity"], m_alfa["Sensitivity"], sens_ci_lo, sens_ci_hi, p_sens, sens_p_boot),
    make_row("Specificity", m_reg["Specificity"], m_alfa["Specificity"], spec_ci_lo, spec_ci_hi, p_spec, spec_p_boot),
    make_row("PPV",         m_reg["PPV"],         m_alfa["PPV"],         ppv_ci_lo,  ppv_ci_hi,  None,   ppv_p_boot),
    make_row("NPV",         m_reg["NPV"],         m_alfa["NPV"],         npv_ci_lo,  npv_ci_hi,  None,   npv_p_boot),
]

res = pd.DataFrame(rows)

# --------------------
# Add BH FDR column (one p-value per metric)
# Use McNemar exact where available; otherwise use bootstrap p
# --------------------
def select_p(row):
    exact = row["Exact p (McNemar)"]
    return exact if pd.notna(exact) else row["Bootstrap p"]

res["Raw p (selected)"] = res.apply(select_p, axis=1).astype(float)

# Benjamini–Hochberg FDR (manual, no external dependency)
# res["FDR_BH"] = bh_fdr(res["Raw p (selected)"].values)

# Optional rounding for readability
for c in ["Exact p (McNemar)", "Bootstrap p", "Raw p (selected)"]:
    res[c] = res[c].astype(float).round(6)

# --------------------
# Save
# --------------------
res.to_csv(OUT_CSV, index=False)

print("\n=== Paired comparison results ===")
print(res.to_string(index=False))
print(f"\nSaved results to: {OUT_CSV}")
