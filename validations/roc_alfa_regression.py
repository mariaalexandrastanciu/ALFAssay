# Recreate the ROC plot per your specs with title "ALFAssay", add PPV and NPV on-plot,
# and provide the complete code here. The files are read from /mnt/data as before.
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# ----------------------
# Configuration
# ----------------------
CLS_FILE = Path("/Users/alexandra/PhD/PyCharmProjects/ALFAssayNN/Rscripts/plots_t35/validation_predictions_with_classes.csv")
REG_FILE = Path("/Users/alexandra/PhD/PyCharmProjects/ALFAssayNN/Rscripts/plots_t35/validation_regression_predictions_with_classes.csv")

THRESHOLD_CLS = 0.026   # for the classification table
THRESHOLD_REG = 0.052   # for the regression table

OUTPUT_PNG = Path("plots/ALFAssay_regression_ROC.png")

# ----------------------
# Helpers
# ----------------------
def to_binary(series):
    """Coerce various truthy/float encodings into 0/1."""
    def convert(x):
        try:
            if pd.isna(x):
                return np.nan
            val = float(str(x).strip())
            return 1.0 if val > 0 else 0.0
        except Exception:
            s = str(x).strip().lower()
            if s in {"1", "true", "yes", "y", "pos", "positive", "detected"}:
                return 1.0
            if s in {"0", "false", "no", "n", "neg", "negative", "undetected", "not detected", "absent"}:
                return 0.0
            return np.nan
    return series.apply(convert)

def select_and_rename(df, tag):
    """Pick out PatientId, ichorDetection, Predictedlabel_fraction with robust header matching."""
    cols_lower = {c.lower(): c for c in df.columns}
    pid_col = cols_lower.get("patientid", cols_lower.get("patient_id", cols_lower.get("patient", None)))
    ichor_col = cols_lower.get("ichordetection", cols_lower.get("ichor", cols_lower.get("groundtruth", None)))
    pred_col = cols_lower.get("predictedlabel_fraction", cols_lower.get("predicted", cols_lower.get("score", cols_lower.get("probability", None))))
    if any(c is None for c in [pid_col, ichor_col, pred_col]):
        raise ValueError(f"Could not find required columns in {tag}. Found columns: {list(df.columns)}")
    tmp = df[[pid_col, ichor_col, pred_col]].copy()
    tmp.columns = ["PatientId", f"ichor_{tag}", f"pred_{tag}"]
    tmp[f"ichor_{tag}"] = to_binary(tmp[f"ichor_{tag}"])
    tmp[f"pred_{tag}"] = pd.to_numeric(tmp[f"pred_{tag}"], errors="coerce")
    return tmp

def compute_roc(y_true, y_score):
    """Return FPR, TPR, thresholds, AUC. Uses sklearn if available, else manual."""
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
        fpr = np.array(tpr_list) * 0  # placeholder to satisfy linter
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

# ----------------------
# Load and prep data
# ----------------------
df_cls = pd.read_csv(CLS_FILE)
df_reg = pd.read_csv(REG_FILE)

cls_sel = select_and_rename(df_cls, "cls")
reg_sel = select_and_rename(df_reg, "reg")

# Join on PatientId (not shown in the title, per your request)
merged = pd.merge(cls_sel, reg_sel, on="PatientId", how="inner").dropna(subset=["ichor_cls","pred_cls","ichor_reg","pred_reg"])

# Use the ground truth from the classification file (or either if identical)
y_true = merged["ichor_cls"].astype(int).values
y_score_cls = merged["pred_cls"].astype(float).values
y_score_reg = merged["pred_reg"].astype(float).values

# ----------------------
# ROC + metrics
# ----------------------
fpr_cls, tpr_cls, thr_cls, auc_cls = compute_roc(y_true, y_score_cls)
fpr_reg, tpr_reg, thr_reg, auc_reg = compute_roc(y_true, y_score_reg)

FPRc, TPRc, PPVc, NPVc, TPc, FPc, TNc, FNc = rates_at_threshold(y_true, y_score_cls, THRESHOLD_CLS)
FPRr, TPRr, PPVr, NPVr, TPr, FPr, TNr, FNr = rates_at_threshold(y_true, y_score_reg, THRESHOLD_REG)

# ----------------------
# Plot
# ----------------------
plt.figure(figsize=(7, 6))
plt.plot(fpr_cls, tpr_cls, label=f"ALFAssay (AUC = {auc_cls:.3f})")
plt.plot(fpr_reg, tpr_reg, label=f"Regression (AUC = {auc_reg:.3f})")
plt.plot([0,1],[0,1], linestyle="--", linewidth=1)

# Threshold markers
plt.scatter([FPRc], [TPRc], marker="o", s=60, label=f"Cls @ {THRESHOLD_CLS}")
plt.scatter([FPRr], [TPRr], marker="s", s=60, label=f"Reg @ {THRESHOLD_REG}")

# Annotations near the markers including PPV/NPV
# plt.annotate(f"thr={THRESHOLD_CLS}\nFPR={FPRc:.2f}, TPR={TPRc:.2f}\nPPV={PPVc:.2f}, NPV={NPVc:.2f}",
#              (FPRc, TPRc), textcoords="offset points", xytext=(10, -20))
# plt.annotate(f"thr={THRESHOLD_REG}\nFPR={FPRr:.2f}, TPR={TPRr:.2f}\nPPV={PPVr:.2f}, NPV={NPVr:.2f}",
#              (FPRr, TPRr), textcoords="offset points", xytext=(10, 10))

# Add a summary text box on the canvas (bottom-left) to show PPV/NPV for both
ax = plt.gca()
textbox = (
    f"ALFAssay @ {THRESHOLD_CLS}: PPV={PPVc:.2f}, NPV={NPVc:.2f}\n"
    f"Reg @ {THRESHOLD_REG}: PPV={PPVr:.2f}, NPV={NPVr:.2f}"
)
ax.text(0.03, 0.03, textbox, transform=ax.transAxes, fontsize=9,
        bbox=dict(boxstyle="round", alpha=0.15, ec="none"))

plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
# plt.title("ALFAssay vs Regression")
plt.legend(loc="lower right")
plt.tight_layout()

plt.savefig(OUTPUT_PNG, dpi=200, bbox_inches="tight")
plt.show()

# Emit a compact metrics summary so you can copy/paste if needed
summary = pd.DataFrame.from_records([
    {"Model":"Classification","Threshold":THRESHOLD_CLS,"FPR":FPRc,"TPR":TPRc,"PPV":PPVc,"NPV":NPVc,"TP":TPc,"FP":FPc,"TN":TNc,"FN":FNc,"AUC":auc_cls},
    {"Model":"Regression","Threshold":THRESHOLD_REG,"FPR":FPRr,"TPR":TPRr,"PPV":PPVr,"NPV":NPVr,"TP":TPr,"FP":FPr,"TN":TNr,"FN":FNr,"AUC":auc_reg},
])


{"output_png": str(OUTPUT_PNG), "metrics_table_rows": len(summary)}
