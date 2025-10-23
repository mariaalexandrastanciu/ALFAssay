# Created by alexandra at 26/09/2025
# UPDATED: Calibrate FRAGLE to ALFAssay scale via linear regression before summing
# 1x3 grid of ROC plots with a 4th curve:
# A: ichorCNA reference (fixed thr for ALFAssay/Reg; Youden for FRAGLE; Youden for ALFA + calibrated FRAGLE sum)
# B: ctDNADetected (NIPT-PRL consistent) reference (Youden for ALL, incl. ALFA + calibrated FRAGLE sum)
# C: PFS reference (pfs_binary: ≤6 mo vs >6 mo) (Youden for ALL, incl. ALFA + calibrated FRAGLE sum)

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# ----------------------
# Configuration
# ----------------------
CLS_FILE    = Path("/Users/alexandra/PhD/PyCharmProjects/ALFAssayNN/Rscripts/plots_t35/validation_predictions_with_classes.csv")
REG_FILE    = Path("/Users/alexandra/PhD/PyCharmProjects/ALFAssayNN/Rscripts/plots_t35/validation_regression_predictions_with_classes.csv")
FRAGLE_FILE = Path("/Users/alexandra/PhD/PyCharmProjects/ALFAssayNN/Rscripts/plots_t35/fragle_predictions_with_classes.csv")
PFS_FILE    = Path("/Users/alexandra/PhD/PyCharmProjects/ALFAssayNN/Rscripts/plots_t35/pfs_pearl.csv")

THRESHOLD_CLS_FIXED = 0.026   # Panel A: fixed ALFAssay threshold
THRESHOLD_REG_FIXED = 0.052   # Panel A: fixed Regression threshold

OUT_GRID = Path("plots/roc_grid.pdf")
OUT_GRID.parent.mkdir(parents=True, exist_ok=True)

# ----------------------
# Helpers
# ----------------------
def to_binary(series):
    def convert(x):
        try:
            if pd.isna(x):
                return np.nan
            val = float(str(x).strip())
            return 1.0 if val > 0 else 0.0
        except Exception:
            s = str(x).strip().lower()
            if s in {"1","true","yes","y","pos","positive","detected"}:
                return 1.0
            if s in {"0","false","no","n","neg","negative","undetected","not detected","absent"}:
                return 0.0
            return np.nan
    return series.apply(convert)

def select_and_rename(df, tag, label_kind):
    cols_lower = {c.lower(): c for c in df.columns}
    pid_col = cols_lower.get("patientid", cols_lower.get("patient_id", cols_lower.get("patient", None)))

    if label_kind == "ichor":
        label_col = cols_lower.get("ichordetection", cols_lower.get("ichor", cols_lower.get("groundtruth", None)))
        label_name = f"ichor_{tag}"
    elif label_kind == "ctdna":
        label_col = cols_lower.get("ctdnadetected", cols_lower.get("ctdna", None))
        label_name = f"ctDNA_{tag}"
    else:
        raise ValueError("label_kind must be 'ichor' or 'ctdna'")

    pred_col = cols_lower.get("predictedlabel_fraction", None)

    if any(c is None for c in [pid_col, label_col, pred_col]):
        raise ValueError(f"Missing required columns in {tag}. Have: {list(df.columns)}")

    tmp = df[[pid_col, label_col, pred_col]].copy()
    tmp.columns = ["PatientId", label_name, f"pred_{tag}"]
    tmp[label_name] = to_binary(tmp[label_name])
    tmp[f"pred_{tag}"] = pd.to_numeric(tmp[f"pred_{tag}"], errors="coerce")
    return tmp

def select_pfs(df):
    cols_lower = {c.lower(): c for c in df.columns}
    pid_col = cols_lower.get("patientid", cols_lower.get("patient_id", cols_lower.get("patient", None)))
    pfs_col = cols_lower.get("pfs_binary", None)
    if any(c is None for c in [pid_col, pfs_col]):
        raise ValueError(f"PFS file must contain patientId and pfs_binary. Found: {list(df.columns)}")
    tmp = df[[pid_col, pfs_col]].copy()
    tmp.columns = ["PatientId", "PFS_binary"]
    tmp["PFS_binary"] = to_binary(tmp["PFS_binary"])
    return tmp

def compute_roc(y_true, y_score):
    try:
        from sklearn.metrics import roc_curve, auc
        fpr, tpr, thr = roc_curve(y_true, y_score)
        roc_auc = auc(fpr, tpr)
        return fpr, tpr, thr, roc_auc
    except Exception:
        thresholds = np.unique(y_score)[::-1]
        thresholds = np.r_[np.inf, thresholds, -np.inf]
        tpr_list, fpr_list = [], []
        for t in thresholds:
            y_pred = (y_score >= t).astype(int)
            TP = ((y_true == 1) & (y_pred == 1)).sum()
            FP = ((y_true == 0) & (y_pred == 1)).sum()
            FN = ((y_true == 1) & (y_pred == 0)).sum()
            TN = ((y_true == 0) & (y_pred == 0)).sum()
            tpr_list.append(TP / (TP + FN) if (TP + FN) > 0 else 0.0)
            fpr_list.append(FP / (FP + TN) if (FP + TN) > 0 else 0.0)
        fpr = np.array(fpr_list)
        tpr = np.array(tpr_list)
        order = np.argsort(fpr)
        roc_auc = np.trapz(tpr[order], fpr[order])
        return fpr, tpr, thresholds, roc_auc

def confusion_at_threshold(y_true, y_score, thr):
    y_pred = (y_score >= thr).astype(int)
    TP = int(((y_true == 1) & (y_pred == 1)).sum())
    FP = int(((y_true == 0) & (y_pred == 1)).sum())
    FN = int(((y_true == 1) & (y_pred == 0)).sum())
    TN = int(((y_true == 0) & (y_pred == 0)).sum())
    return TP, FP, TN, FN

def rates_at_threshold(y_true, y_score, thr):
    TP, FP, TN, FN = confusion_at_threshold(y_true, y_score, thr)
    TPR = TP / (TP + FN) if (TP + FN) else np.nan
    FPR = FP / (FP + TN) if (FP + TN) else np.nan
    PPV = TP / (TP + FP) if (TP + FP) else np.nan
    NPV = TN / (TN + FN) if (TN + FN) else np.nan
    return FPR, TPR, PPV, NPV, TP, FP, TN, FN

def youden_cutoff(fpr, tpr, thr):
    J = tpr - fpr
    idx = int(np.nanargmax(J))
    return float(thr[idx]), float(J[idx]), idx

def add_panel_label(ax, label):
    ax.text(-0.1, 1.05, f"{label}.", transform=ax.transAxes,
            ha="left", va="bottom", fontsize=12, fontweight="bold")

# --- NEW: FRAGLE -> ALFA calibration ---
def fit_linear_calibration(x, y):
    """
    Fit y ≈ a + b*x. Returns (a, b).
    Tries sklearn LinearRegression, falls back to numpy polyfit.
    """
    x = np.asarray(x).reshape(-1, 1)
    y = np.asarray(y).reshape(-1, 1)
    try:
        from sklearn.linear_model import LinearRegression
        mdl = LinearRegression()  # with intercept
        mdl.fit(x, y)
        a = float(mdl.intercept_[0])
        b = float(mdl.coef_[0, 0])
    except Exception:
        b, a = np.polyfit(x.ravel(), y.ravel(), 1)  # returns slope, intercept
    return a, b

def apply_calibration(x, a, b, clip=True):
    y = a + b * np.asarray(x, dtype=float)
    if clip:
        y = np.clip(y, 0.0, 1.0)
    return y

# ----------------------
# Load data once
# ----------------------
df_cls    = pd.read_csv(CLS_FILE)     # ALFAssay classification file (Predictedlabel_fraction)
df_reg    = pd.read_csv(REG_FILE)     # Regression file
df_fragle = pd.read_csv(FRAGLE_FILE)  # FRAGLE file (Predictedlabel_fraction)
df_pfs    = pd.read_csv(PFS_FILE)

# ---------- Panel A: ichorCNA reference ----------
cls_ichor = select_and_rename(df_cls, "cls", "ichor")
reg_ichor = select_and_rename(df_reg, "reg", "ichor")
fra_ichor = select_and_rename(df_fragle, "fragle", "ichor")

merged_A = (
    cls_ichor.merge(reg_ichor, on="PatientId", how="inner")
             .merge(fra_ichor,  on="PatientId", how="inner")
             .dropna(subset=["ichor_cls","pred_cls","ichor_reg","pred_reg","pred_fragle"])
)

# --- NEW: Fit calibration on ALFA vs FRAGLE predictions (global), using all patients with both ---
calib_df = (
    cls_ichor[["PatientId", "pred_cls"]]
    .merge(fra_ichor[["PatientId", "pred_fragle"]], on="PatientId", how="inner")
    .dropna(subset=["pred_cls", "pred_fragle"])
)
if calib_df.empty:
    raise ValueError("No overlap between ALFAssay and FRAGLE predictions to fit calibration.")
a_cal, b_cal = fit_linear_calibration(calib_df["pred_fragle"].values, calib_df["pred_cls"].values)
print(f"[Calibration] FRAGLE_to_ALFA: ALFA ≈ a + b*FRAGLE with a={a_cal:.4f}, b={b_cal:.4f}")

# Use ichor_cls as ground truth for A
y_true_A = merged_A["ichor_cls"].astype(int).values
y_cls_A  = merged_A["pred_cls"].astype(float).values
y_reg_A  = merged_A["pred_reg"].astype(float).values
y_fra_A  = merged_A["pred_fragle"].astype(float).values

# --- NEW: Calibrated FRAGLE and calibrated sum ---
y_fra_cal_A = apply_calibration(y_fra_A, a_cal, b_cal, clip=True)
y_sum_A     = y_cls_A + y_fra_cal_A

# ROCs
fpr_cls_A, tpr_cls_A, thr_cls_A, auc_cls_A = compute_roc(y_true_A, y_cls_A)
fpr_reg_A, tpr_reg_A, thr_reg_A, auc_reg_A = compute_roc(y_true_A, y_reg_A)
fpr_fra_A, tpr_fra_A, thr_fra_A, auc_fra_A = compute_roc(y_true_A, y_fra_A)
fpr_sum_A, tpr_sum_A, thr_sum_A, auc_sum_A = compute_roc(y_true_A, y_sum_A)

# Thresholds (Panel A): fixed for ALFA/Reg; Youden for FRAGLE and SUM
FPRc_A, TPRc_A, PPVc_A, NPVc_A, TPc_A, FPc_A, TNc_A, FNc_A = rates_at_threshold(y_true_A, y_cls_A, THRESHOLD_CLS_FIXED)
FPRr_A, TPRr_A, PPVr_A, NPVr_A, TPr_A, FPr_A, TNr_A, FNr_A = rates_at_threshold(y_true_A, y_reg_A, THRESHOLD_REG_FIXED)
thr_fra_opt_A, J_fra_A, _ = youden_cutoff(fpr_fra_A, tpr_fra_A, thr_fra_A)
FPRf_A, TPRf_A, PPVf_A, NPVf_A, TPf_A, FPf_A, TNf_A, FNf_A = rates_at_threshold(y_true_A, y_fra_A, thr_fra_opt_A)
thr_sum_opt_A, J_sum_A, _ = youden_cutoff(fpr_sum_A, tpr_sum_A, thr_sum_A)
FPRs_A, TPRs_A, PPVs_A, NPVs_A, TPs_A, FPs_A, TNs_A, FNs_A = rates_at_threshold(y_true_A, y_sum_A, thr_sum_opt_A)

# ---------- Panel B: ctDNA reference (NIPT-PRL only; all equal) ----------
cls_ct = select_and_rename(df_cls, "cls", "ctdna")
reg_ct = select_and_rename(df_reg, "reg", "ctdna")
fra_ct = select_and_rename(df_fragle, "fragle", "ctdna")

merged_B = (
    cls_ct.merge(reg_ct, on="PatientId", how="inner")
          .merge(fra_ct,  on="PatientId", how="inner")
)
merged_B = merged_B.dropna(subset=["ctDNA_cls","ctDNA_reg","ctDNA_fragle"])
mask = merged_B["PatientId"].astype(str).str.startswith("NIPT-PRL", na=False)
merged_B = merged_B[mask].copy()

merged_B["all_equal"] = ((merged_B["ctDNA_cls"] == merged_B["ctDNA_reg"]) &
                         (merged_B["ctDNA_cls"] == merged_B["ctDNA_fragle"])).astype(int)
data_B = merged_B.loc[merged_B["all_equal"] == 1].copy().dropna(subset=["pred_cls","pred_reg","pred_fragle"])

y_true_B = data_B["ctDNA_cls"].astype(int).values
y_cls_B  = data_B["pred_cls"].astype(float).values
y_reg_B  = data_B["pred_reg"].astype(float).values
y_fra_B  = data_B["pred_fragle"].astype(float).values

# --- NEW: Calibrated FRAGLE and calibrated sum (Panel B) ---
y_fra_cal_B = apply_calibration(y_fra_B, a_cal, b_cal, clip=True)
y_sum_B     = y_cls_B + y_fra_cal_B

# ROCs
fpr_cls_B, tpr_cls_B, thr_cls_B, auc_cls_B = compute_roc(y_true_B, y_cls_B)
fpr_reg_B, tpr_reg_B, thr_reg_B, auc_reg_B = compute_roc(y_true_B, y_reg_B)
fpr_fra_B, tpr_fra_B, thr_fra_B, auc_fra_B = compute_roc(y_true_B, y_fra_B)
fpr_sum_B, tpr_sum_B, thr_sum_B, auc_sum_B = compute_roc(y_true_B, y_sum_B)

# Thresholds = Youden for all in Panel B
thr_cls_opt_B, J_cls_B, _ = youden_cutoff(fpr_cls_B, tpr_cls_B, thr_cls_B)
thr_reg_opt_B, J_reg_B, _ = youden_cutoff(fpr_reg_B, tpr_reg_B, thr_reg_B)
thr_fra_opt_B, J_fra_B, _ = youden_cutoff(fpr_fra_B, tpr_fra_B, thr_fra_B)
thr_sum_opt_B, J_sum_B, _ = youden_cutoff(fpr_sum_B, tpr_sum_B, thr_sum_B)

FPRc_B, TPRc_B, PPVc_B, NPVc_B, TPc_B, FPc_B, TNc_B, FNc_B = rates_at_threshold(y_true_B, y_cls_B, thr_cls_opt_B)
FPRr_B, TPRr_B, PPVr_B, NPVr_B, TPr_B, FPr_B, TNr_B, FNr_B = rates_at_threshold(y_true_B, y_reg_B, thr_reg_opt_B)
FPRf_B, TPRf_B, PPVf_B, NPVf_B, TPf_B, FPf_B, TNf_B, FNf_B = rates_at_threshold(y_true_B, y_fra_B, thr_fra_opt_B)
FPRs_B, TPRs_B, PPVs_B, NPVs_B, TPs_B, FPs_B, TNs_B, FNs_B = rates_at_threshold(y_true_B, y_sum_B, thr_sum_opt_B)

# ---------- Panel C: PFS reference ----------
pfs_sel = select_pfs(df_pfs)

preds_all = (
    cls_ichor[["PatientId","pred_cls"]]
      .merge(reg_ichor[["PatientId","pred_reg"]], on="PatientId", how="inner")
      .merge(fra_ichor[["PatientId","pred_fragle"]], on="PatientId", how="inner")
)
merged_C = preds_all.merge(pfs_sel, on="PatientId", how="inner").dropna(subset=["PFS_binary","pred_cls","pred_reg","pred_fragle"])
if merged_C.empty:
    raise ValueError("No rows available for PFS ROC (check joins and non-null PFS_binary).")

y_true_C = merged_C["PFS_binary"].astype(int).values
y_cls_C  = merged_C["pred_cls"].astype(float).values
y_reg_C  = merged_C["pred_reg"].astype(float).values
y_fra_C  = merged_C["pred_fragle"].astype(float).values

# --- NEW: Calibrated FRAGLE and calibrated sum (Panel C) ---
y_fra_cal_C = apply_calibration(y_fra_C, a_cal, b_cal, clip=True)
y_sum_C     = y_cls_C + y_fra_cal_C

# ROCs
fpr_cls_C, tpr_cls_C, thr_cls_C, auc_cls_C = compute_roc(y_true_C, y_cls_C)
fpr_reg_C, tpr_reg_C, thr_reg_C, auc_reg_C = compute_roc(y_true_C, y_reg_C)
fpr_fra_C, tpr_fra_C, thr_fra_C, auc_fra_C = compute_roc(y_true_C, y_fra_C)
fpr_sum_C, tpr_sum_C, thr_sum_C, auc_sum_C = compute_roc(y_true_C, y_sum_C)

# Thresholds = Youden for all in Panel C
thr_cls_opt_C, J_cls_C, _ = youden_cutoff(fpr_cls_C, tpr_cls_C, thr_cls_C)
thr_reg_opt_C, J_reg_C, _ = youden_cutoff(fpr_reg_C, tpr_reg_C, thr_reg_C)
thr_fra_opt_C, J_fra_C, _ = youden_cutoff(fpr_fra_C, tpr_fra_C, thr_fra_C)
thr_sum_opt_C, J_sum_C, _ = youden_cutoff(fpr_sum_C, tpr_sum_C, thr_sum_C)

FPRc_C, TPRc_C, PPVc_C, NPVc_C, TPc_C, FPc_C, TNc_C, FNc_C = rates_at_threshold(y_true_C, y_cls_C, thr_cls_opt_C)
FPRr_C, TPRr_C, PPVr_C, NPVr_C, TPr_C, FPr_C, TNr_C, FNr_C = rates_at_threshold(y_true_C, y_reg_C, thr_reg_opt_C)
FPRf_C, TPRf_C, PPVf_C, NPVf_C, TPf_C, FPf_C, TNf_C, FNf_C = rates_at_threshold(y_true_C, y_fra_C, thr_fra_opt_C)
FPRs_C, TPRs_C, PPVs_C, NPVs_C, TPs_C, FPs_C, TNs_C, FNs_C = rates_at_threshold(y_true_C, y_sum_C, thr_sum_opt_C)

# ----------------------
# Plot all three panels
# ----------------------
fig, axes = plt.subplots(1, 3, figsize=(19, 6), constrained_layout=True)

# ---- Panel A ----
ax = axes[0]
ax.plot(fpr_cls_A, tpr_cls_A, label=f"ALFAssay (AUC = {auc_cls_A:.3f})")
ax.plot(fpr_reg_A, tpr_reg_A, label=f"Regression (AUC = {auc_reg_A:.3f})")
ax.plot(fpr_fra_A, tpr_fra_A, label=f"FRAGLE (AUC = {auc_fra_A:.3f})")
ax.plot(fpr_sum_A, tpr_sum_A, label=f"ALFAssay + FRAGLE (AUC = {auc_sum_A:.3f})")
ax.plot([0,1], [0,1], linestyle="--", linewidth=1, color="gray", alpha=0.5)

ax.scatter([FPRc_A], [TPRc_A], marker="o", s=60, label=f"ALFAssay @ {THRESHOLD_CLS_FIXED}")
ax.scatter([FPRr_A], [TPRr_A], marker="s", s=60, label=f"Reg @ {THRESHOLD_REG_FIXED}")
ax.scatter([FPRf_A], [TPRf_A], marker="^", s=70, label=f"FRAGLE @ {thr_fra_opt_A:.3f} ")
ax.scatter([FPRs_A], [TPRs_A], marker="D", s=70, label=f"ALFAssay+FRAGLE @ {thr_sum_opt_A:.3f}")

ax.set_title("A.", loc="left", fontweight="bold")
ax.set_xlabel("False Positive Rate")
ax.set_ylabel("True Positive Rate")
txtA = (f"ALFAssay @ {THRESHOLD_CLS_FIXED}: PPV={PPVc_A:.2f}, NPV={NPVc_A:.2f}\n"
        f"Reg @ {THRESHOLD_REG_FIXED}: PPV={PPVr_A:.2f}, NPV={NPVr_A:.2f}\n"
        f"FRAGLE @ {thr_fra_opt_A:.3f}: PPV={PPVf_A:.2f}, NPV={NPVf_A:.2f}\n"
        f"ALFAssay+FRAGLE @ {thr_sum_opt_A:.3f}: PPV={PPVs_A:.2f}, NPV={NPVs_A:.2f}")
ax.text(0.03, 0.03, txtA, transform=ax.transAxes, fontsize=9,
        bbox=dict(boxstyle="round", alpha=0.15, ec="none"))
ax.legend(loc="lower right", fontsize=9,
          title="Reference: ichorCNA ctDNA detection", title_fontsize=9)

# ---- Panel B ----
ax = axes[1]
ax.plot(fpr_cls_B, tpr_cls_B, label=f"ALFAssay (AUC = {auc_cls_B:.3f})")
ax.plot(fpr_reg_B, tpr_reg_B, label=f"Regression (AUC = {auc_reg_B:.3f})")
ax.plot(fpr_fra_B, tpr_fra_B, label=f"FRAGLE (AUC = {auc_fra_B:.3f})")
ax.plot(fpr_sum_B, tpr_sum_B, label=f"ALFAssay + FRAGLE (AUC = {auc_sum_B:.3f})")
ax.plot([0,1], [0,1], linestyle="--", linewidth=1, color="gray", alpha=0.5)

ax.scatter([FPRc_B], [TPRc_B], marker="o", s=60, label=f"ALFAssay @ {thr_cls_opt_B:.3f}")
ax.scatter([FPRr_B], [TPRr_B], marker="s", s=60, label=f"Reg @ {thr_reg_opt_B:.3f} ")
ax.scatter([FPRf_B], [TPRf_B], marker="^", s=70, label=f"FRAGLE @ {thr_fra_opt_B:.3f}")
ax.scatter([FPRs_B], [TPRs_B], marker="D", s=70, label=f"ALFAssay+FRAGLE @ {thr_sum_opt_B:.3f} ")

ax.set_title("B.", loc="left", fontweight="bold")
ax.set_xlabel("False Positive Rate")
ax.set_ylabel("True Positive Rate")
txtB = (f"ALFAssay @ {thr_cls_opt_B:.3f}: PPV={PPVc_B:.2f}, NPV={NPVc_B:.2f}\n"
        f"Reg @ {thr_reg_opt_B:.3f}: PPV={PPVr_B:.2f}, NPV={NPVr_B:.2f}\n"
        f"FRAGLE @ {thr_fra_opt_B:.3f}: PPV={PPVf_B:.2f}, NPV={NPVf_B:.2f}\n"
        f"ALFAssay+FRAGLE @ {thr_sum_opt_B:.3f}: PPV={PPVs_B:.2f}, NPV={NPVs_B:.2f}")
ax.text(0.03, 0.03, txtB, transform=ax.transAxes, fontsize=9,
        bbox=dict(boxstyle="round", alpha=0.15, ec="none"))
ax.legend(loc="lower right", fontsize=9,
          title="Reference: OncoFollowDNA ctDNA detection", title_fontsize=9)

# ---- Panel C ----
ax = axes[2]
ax.plot(fpr_cls_C, tpr_cls_C, label=f"ALFAssay (AUC = {auc_cls_C:.3f})")
ax.plot(fpr_reg_C, tpr_reg_C, label=f"Regression (AUC = {auc_reg_C:.3f})")
ax.plot(fpr_fra_C, tpr_fra_C, label=f"FRAGLE (AUC = {auc_fra_C:.3f})")
ax.plot(fpr_sum_C, tpr_sum_C, label=f"ALFAssay + FRAGLE (AUC = {auc_sum_C:.3f})")
ax.plot([0,1], [0,1], linestyle="--", linewidth=1, color="gray", alpha=0.5)

ax.scatter([FPRc_C], [TPRc_C], marker="o", s=60, label=f"ALFAssay @ {thr_cls_opt_C:.3f} ")
ax.scatter([FPRr_C], [TPRr_C], marker="s", s=60, label=f"Reg @ {thr_reg_opt_C:.3f}")
ax.scatter([FPRf_C], [TPRf_C], marker="^", s=70, label=f"FRAGLE @ {thr_fra_opt_C:.3f}")
ax.scatter([FPRs_C], [TPRs_C], marker="D", s=70, label=f"ALFAssay+FRAGLE @ {thr_sum_opt_C:.3f}")

ax.set_title("C.", loc="left", fontweight="bold")
ax.set_xlabel("False Positive Rate")
ax.set_ylabel("True Positive Rate")
txtC = ( f"ALFAssay @ {thr_cls_opt_C:.3f}: PPV={PPVc_C:.2f}, NPV={NPVc_C:.2f}\n"
        f"Reg @ {thr_reg_opt_C:.3f}: PPV={PPVr_C:.2f}, NPV={NPVr_C:.2f}\n"
        f"FRAGLE @ {thr_fra_opt_C:.3f}: PPV={PPVf_C:.2f}, NPV={NPVf_C:.2f}\n"
        f"ALFAssay+FRAGLE @ {thr_sum_opt_C:.3f}: PPV={PPVs_C:.2f}, NPV={NPVs_C:.2f}" )
ax.text(0.03, 0.03, txtC, transform=ax.transAxes, fontsize=9,
        bbox=dict(boxstyle="round", alpha=0.15, ec="none"))
ax.legend(loc="lower right", fontsize=9,
          title="Reference: PFS ≤ 6 mo vs > 6 mo", title_fontsize=9)

# Save grid
fig.savefig(OUT_GRID, dpi=200, bbox_inches="tight")
print(f"Saved grid to: {OUT_GRID.resolve()}")
