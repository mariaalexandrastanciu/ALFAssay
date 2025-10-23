# ROC for ALFAssay, Regression, and FRAGLE using ctDNADetected as ground truth,
# restricted to NIPT-PRL patients with non-null ctDNADetected in all files.
# Thresholds for ALL models (ALFAssay, Regression, FRAGLE) are recalculated via Youden's J.

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# ----------------------
# Configuration
# ----------------------
CLS_FILE = Path("/Users/alexandra/PhD/PyCharmProjects/ALFAssayNN/Rscripts/plots_t35/validation_predictions_with_classes.csv")
REG_FILE = Path("/Users/alexandra/PhD/PyCharmProjects/ALFAssayNN/Rscripts/plots_t35/validation_regression_predictions_with_classes.csv")
FRAGLE_FILE = Path("/Users/alexandra/PhD/PyCharmProjects/ALFAssayNN/Rscripts/plots_t35/fragle_predictions_with_classes.csv")

OUTPUT_PNG = Path("plots/roc_oncodna.png")
INCONSISTENCIES_CSV = Path("plots/ctDNA_inconsistencies_NIPT_PRL.csv")

# ----------------------
# Helpers
# ----------------------
def to_binary(series):
    """Coerce truthy/float encodings into 0/1 (NaN preserved)."""
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

def select_and_rename(df, tag):
    """
    Extract PatientId, ctDNADetected (ctDNA), Predictedlabel_fraction (score).
    Handles common header variants.
    """
    cols_lower = {c.lower(): c for c in df.columns}
    pid_col = cols_lower.get("patientid", cols_lower.get("patient_id", cols_lower.get("patient", None)))
    ctdna_col = cols_lower.get("ctdnadetected", cols_lower.get("ctdna", None))
    pred_col = cols_lower.get("predictedlabel_fraction", cols_lower.get("predicted", cols_lower.get("score", cols_lower.get("probability", None))))
    if any(c is None for c in [pid_col, ctdna_col, pred_col]):
        raise ValueError(f"Missing columns in {tag}. Found: {list(df.columns)}")
    tmp = df[[pid_col, ctdna_col, pred_col]].copy()
    tmp.columns = ["PatientId", f"ctDNA_{tag}", f"pred_{tag}"]
    tmp[f"ctDNA_{tag}"] = to_binary(tmp[f"ctDNA_{tag}"])
    tmp[f"pred_{tag}"] = pd.to_numeric(tmp[f"pred_{tag}"], errors="coerce")
    return tmp

def compute_roc(y_true, y_score):
    from sklearn.metrics import roc_curve, auc
    fpr, tpr, thr = roc_curve(y_true, y_score)
    roc_auc = auc(fpr, tpr)
    return fpr, tpr, thr, roc_auc

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
    """Return threshold that maximizes Youden's J = TPR - FPR."""
    J = tpr - fpr
    idx = int(np.nanargmax(J))
    return float(thr[idx]), float(J[idx]), int(idx)

# ----------------------
# Load and prep data
# ----------------------
df_cls = pd.read_csv(CLS_FILE)
df_reg = pd.read_csv(REG_FILE)
df_fragle = pd.read_csv(FRAGLE_FILE)

cls_sel = select_and_rename(df_cls, "cls")
reg_sel = select_and_rename(df_reg, "reg")
fragle_sel = select_and_rename(df_fragle, "fragle")

# Merge on PatientId
merged = (
    cls_sel
    .merge(reg_sel, on="PatientId", how="inner")
    .merge(fragle_sel, on="PatientId", how="inner")
)

# Keep only rows with non-null ctDNA in ALL files and PatientId starting with NIPT-PRL
merged = merged.dropna(subset=["ctDNA_cls","ctDNA_reg","ctDNA_fragle"])
merged = merged[merged["PatientId"].astype(str).str.startswith("NIPT-PRL")]

# ----------------------
# Consistency check
# ----------------------
merged["all_equal"] = (
    (merged["ctDNA_cls"] == merged["ctDNA_reg"]) &
    (merged["ctDNA_cls"] == merged["ctDNA_fragle"])
).astype(int)

n_total = len(merged)
n_all_equal = int(merged["all_equal"].sum())
print("=== ctDNADetected consistency (NIPT-PRL only) ===")
print(f"Total patients: {n_total}")
print(f"All equal: {n_all_equal} ({(n_all_equal / n_total if n_total else 0):.1%})")

# Save disagreements for inspection
disagreements = merged.loc[merged["all_equal"] == 0, ["PatientId", "ctDNA_cls","ctDNA_reg","ctDNA_fragle"]]
if not disagreements.empty:
    INCONSISTENCIES_CSV.parent.mkdir(parents=True, exist_ok=True)
    disagreements.to_csv(INCONSISTENCIES_CSV, index=False)
    print(f"Disagreements saved to {INCONSISTENCIES_CSV} (rows={len(disagreements)})")
else:
    print("No disagreements across the three files.")

# Use only consistent rows for ROC and drop null predictions
data = merged.loc[merged["all_equal"] == 1].copy().dropna(subset=["pred_cls","pred_reg","pred_fragle"])
print(f"Using {len(data)} consistent NIPT-PRL patients for ROC.")

# ----------------------
# Ground truth + predictions
# ----------------------
y_true = data["ctDNA_cls"].astype(int).values
y_score_cls = data["pred_cls"].astype(float).values
y_score_reg = data["pred_reg"].astype(float).values
y_score_fragle = data["pred_fragle"].astype(float).values

# ----------------------
# ROC + AUC
# ----------------------
fpr_cls, tpr_cls, thr_cls, auc_cls = compute_roc(y_true, y_score_cls)
fpr_reg, tpr_reg, thr_reg, auc_reg = compute_roc(y_true, y_score_reg)
fpr_fra, tpr_fra, thr_fra, auc_fra = compute_roc(y_true, y_score_fragle)

# ----------------------
# Recalculate thresholds via Youden's J for ALL models
# ----------------------
thr_cls_opt, J_cls, _ = youden_cutoff(fpr_cls, tpr_cls, thr_cls)
thr_reg_opt, J_reg, _ = youden_cutoff(fpr_reg, tpr_reg, thr_reg)
thr_fra_opt, J_fra, _ = youden_cutoff(fpr_fra, tpr_fra, thr_fra)

# Metrics at optimal thresholds
FPRc, TPRc, PPVc, NPVc, TPc, FPc, TNc, FNc = rates_at_threshold(y_true, y_score_cls, thr_cls_opt)
FPRr, TPRr, PPVr, NPVr, TPr, FPr, TNr, FNr = rates_at_threshold(y_true, y_score_reg, thr_reg_opt)
FPRf, TPRf, PPVf, NPVf, TPf, FPf, TNf, FNf = rates_at_threshold(y_true, y_score_fragle, thr_fra_opt)

# ----------------------
# Plot
# ----------------------
plt.figure(figsize=(7, 6))
plt.plot(fpr_cls, tpr_cls, label=f"ALFAssay (AUC = {auc_cls:.3f})")
plt.plot(fpr_reg, tpr_reg, label=f"Regression (AUC = {auc_reg:.3f})")
plt.plot(fpr_fra, tpr_fra, label=f"FRAGLE (AUC = {auc_fra:.3f})")
plt.plot([0,1],[0,1], linestyle="--", linewidth=1)

# Threshold markers at Youden-optimal thresholds for ALL models
plt.scatter([FPRc], [TPRc], marker="o", s=60, label=f"ALFAssay @ {thr_cls_opt:.3f} (Youden)")
plt.scatter([FPRr], [TPRr], marker="s", s=60, label=f"Reg @ {thr_reg_opt:.3f} (Youden)")
plt.scatter([FPRf], [TPRf], marker="^", s=70, label=f"FRAGLE @ {thr_fra_opt:.3f} (Youden)")

# Info textbox with PPV/NPV for all three at the optimal thresholds
ax = plt.gca()
textbox = (
    f"ALFAssay @ {thr_cls_opt:.3f}: PPV={PPVc:.2f}, NPV={NPVc:.2f}\n"
    f"Regression @ {thr_reg_opt:.3f}: PPV={PPVr:.2f}, NPV={NPVr:.2f}\n"
    f"FRAGLE @ {thr_fra_opt:.3f}: PPV={PPVf:.2f}, NPV={NPVf:.2f}"
)
ax.text(0.03, 0.03, textbox, transform=ax.transAxes, fontsize=9,
        bbox=dict(boxstyle="round", alpha=0.15, ec="none"))

plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("ROC on oncoDNA as orthogonal validation")
plt.legend(loc="lower right")
plt.tight_layout()

OUTPUT_PNG.parent.mkdir(parents=True, exist_ok=True)
plt.savefig(OUTPUT_PNG, dpi=200, bbox_inches="tight")
plt.show()

# ----------------------
# Summary table
# ----------------------
summary = pd.DataFrame.from_records([
    {"Model":"ALFAssay (Cls)","Threshold_Youden":thr_cls_opt,"FPR":FPRc,"TPR":TPRc,"PPV":PPVc,"NPV":NPVc,"TP":TPc,"FP":FPc,"TN":TNc,"FN":FNc,"AUC":auc_cls,"Youden_J":J_cls},
    {"Model":"Regression","Threshold_Youden":thr_reg_opt,"FPR":FPRr,"TPR":TPRr,"PPV":PPVr,"NPV":NPVr,"TP":TPr,"FP":FPr,"TN":TNr,"FN":FNr,"AUC":auc_reg,"Youden_J":J_reg},
    {"Model":"FRAGLE","Threshold_Youden":thr_fra_opt,"FPR":FPRf,"TPR":TPRf,"PPV":PPVf,"NPV":NPVf,"TP":TPf,"FP":FPf,"TN":TNf,"FN":FNf,"AUC":auc_fra,"Youden_J":J_fra},
])

print("\n=== ROC Summary (NIPT-PRL, Youden thresholds) ===")
print(summary.to_string(index=False))
