# Created by alexandra at 12/08/2025
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc, confusion_matrix, accuracy_score, recall_score, precision_score

# --------------------
# Config
# --------------------
TRAIN_PATH ="/Users/alexandra/PhD/PyCharmProjects/ALFAssayNN/Rscripts/plots_t35/data_regression_training_ctDNADetected.csv"
         # TSV
VAL_PATH   = "/Users/alexandra/PhD/PyCharmProjects/ALFAssayNN/Rscripts/plots_t35/data_regression_validation_ctDNADetected.csv" # TSV
N_BOOT = 2000
RANDOM_STATE = 42
OUT_TRAIN_PRED = "/Users/alexandra/PhD/PyCharmProjects/ALFAssayNN/Rscripts/plots_t35/training_regression_predictions_with_classes.csv"
OUT_VAL_PRED   = "/Users/alexandra/PhD/PyCharmProjects/ALFAssayNN/Rscripts/plots_t35/validation_regression_predictions_with_classes.csv"
OUT_TRAIN_ROC  = "plots/cut_off_regression_training_roc.png"
OUT_VAL_ROC  = "plots/cut_off_regression_validation_roc.png"
OUT_METRICS    =  "/Users/alexandra/PhD/PyCharmProjects/ALFAssayNN/Rscripts/plots_t35/regressionModel_metrics_summary.csv"
rng = np.random.default_rng(RANDOM_STATE)

# --------------------
# Load + clean
# --------------------
train = pd.read_csv(TRAIN_PATH, sep="\t")
val   = pd.read_csv(VAL_PATH,   sep="\t")

# Remove rows with null ctDNADetected
train = train.dropna(subset=["ichorDetection"]).copy()
val   = val.dropna(subset=["ichorDetection"]).copy()

# Convert Predictedlabel (percentage) -> fraction [0,1]
train["Predictedlabel_fraction"] = train["prediction"]
val["Predictedlabel_fraction"]   = val["prediction"]

# Ensure ground truth is binary int {0,1}
train["ichorDetection"] = train["ichorDetection"].astype(int)
val["ichorDetection"]   = val["ichorDetection"].astype(int)

# --------------------
# Helper functions
# --------------------
def youden_cutoff(y_true, y_score):
    fpr, tpr, thr = roc_curve(y_true, y_score)
    j = tpr - fpr
    idx = np.argmax(j)
    return thr[idx], (fpr, tpr, thr), auc(fpr, tpr), idx

def metrics_at_cutoff(y_true, y_score, cutoff):
    """Return confusion counts + metrics (Acc, Sens, Spec, PPV, NPV) at cutoff."""
    y_pred = (y_score >= cutoff).astype(int)
    # Force labels=[0,1] so we always get 2x2, even if a class is missing in a bootstrap sample
    cm = confusion_matrix(y_true, y_pred, labels=[0, 1])
    tn, fp, fn, tp = cm.ravel()
    acc = accuracy_score(y_true, y_pred)
    sens = recall_score(y_true, y_pred)  # TPR
    spec = recall_score(y_true, y_pred, pos_label=0)  # TNR
    ppv = tp / (tp + fp) if (tp + fp) > 0 else np.nan
    npv = tn / (tn + fn) if (tn + fn) > 0 else np.nan
    # precision == PPV; we compute PPV explicitly above and keep precision for completeness
    prec = precision_score(y_true, y_pred, zero_division=0)
    return {
        "tn": tn, "fp": fp, "fn": fn, "tp": tp,
        "accuracy": acc, "sensitivity": sens, "specificity": spec,
        "precision": prec, "ppv": ppv, "npv": npv, "cm": cm
    }


def ci_from_boot(arr, alpha=0.05):
    """Percentile bootstrap CI."""
    lo = np.nanpercentile(arr, 100 * alpha / 2)
    hi = np.nanpercentile(arr, 100 * (1 - alpha / 2))
    return float(lo), float(hi)

# --------------------
# Train: choose cutoff (Youden)
# --------------------
y_true_tr = train["ichorDetection"].to_numpy()
y_score_tr = train["Predictedlabel_fraction"].to_numpy()

cutoff, (fpr_tr, tpr_tr, thr_tr), auc_tr, idx_tr = youden_cutoff(y_true_tr, y_score_tr)
print(f"Training AUC: {auc_tr:.3f}")
print(f"Chosen cut-off (fraction): {cutoff:.6f}  |  ({cutoff*100:.3f}%)")

# Plot training ROC
plt.figure()
plt.plot(fpr_tr, tpr_tr, label=f"Train ROC (AUC = {auc_tr:.3f})")
plt.scatter(fpr_tr[idx_tr], tpr_tr[idx_tr], label=f"Cut-off = {cutoff:.4f}", zorder=5)
plt.plot([0, 1], [0, 1], 'k--')
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("Training ROC: Predictedlabel vs ctDNADetected")
plt.legend()
plt.tight_layout()
plt.savefig(OUT_TRAIN_ROC, dpi=150)
plt.close()

# Training metrics at chosen cutoff (useful for the summary table)
m_tr = metrics_at_cutoff(y_true_tr, y_score_tr, cutoff)

# --------------------
# Validation: metrics with the same cutoff
# --------------------
y_true_val = val["ichorDetection"].to_numpy()
y_score_val = val["Predictedlabel_fraction"].to_numpy()

m_val = metrics_at_cutoff(y_true_val, y_score_val, cutoff)
print("\nValidation performance (using training cut-off):")
print("Confusion matrix [[TN, FP], [FN, TP]]:")
print(m_val["cm"])
print(f"Accuracy   : {m_val['accuracy']:.4f}")
print(f"Sensitivity: {m_val['sensitivity']:.4f}")
print(f"Specificity: {m_val['specificity']:.4f}")
print(f"PPV        : {m_val['ppv']:.4f}")
print(f"NPV        : {m_val['npv']:.4f}")

# Validation ROC (mark the chosen cutoff point)
fpr_v, tpr_v, thr_v = roc_curve(y_true_val, y_score_val)
auc_v = auc(fpr_v, tpr_v)
closest_idx = np.argmin(np.abs(thr_v - cutoff))

plt.figure()
plt.plot(fpr_v, tpr_v, label=f"Validation ROC (AUC = {auc_v:.3f})")
plt.scatter(fpr_v[closest_idx], tpr_v[closest_idx], label=f"Cut-off = {cutoff:.4f}", zorder=5)
plt.plot([0, 1], [0, 1], 'k--')
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("Validation ROC: Predictedlabel vs ctDNADetected")
plt.legend()
plt.tight_layout()
plt.savefig(OUT_VAL_ROC, dpi=150)
plt.close()

# --------------------
# Bootstrap CI for cutoff (training)
# --------------------
cut_boot = np.empty(N_BOOT, dtype=float)
n_tr = len(y_true_tr)

for b in range(N_BOOT):
    idx = rng.integers(0, n_tr, size=n_tr)   # resample with replacement
    yb = y_true_tr[idx]
    sb = y_score_tr[idx]
    try:
        c_b, _, _, _ = youden_cutoff(yb, sb)
    except Exception:
        c_b = np.nan
    cut_boot[b] = c_b

cut_boot = cut_boot[~np.isnan(cut_boot)]
cut_ci = ci_from_boot(cut_boot, alpha=0.05)
print(f"\nBootstrap 95% CI for cut-off (fraction): [{cut_ci[0]:.6f}, {cut_ci[1]:.6f}]")
print(f"Bootstrap 95% CI for cut-off (percent) : [{cut_ci[0]*100:.3f}%, {cut_ci[1]*100:.3f}%]")

# --------------------
# Bootstrap CIs for validation metrics (fixed cutoff)
# --------------------
n_val = len(y_true_val)
acc_b  = np.empty(N_BOOT, dtype=float)
sens_b = np.empty(N_BOOT, dtype=float)
spec_b = np.empty(N_BOOT, dtype=float)
ppv_b  = np.empty(N_BOOT, dtype=float)
npv_b  = np.empty(N_BOOT, dtype=float)

for b in range(N_BOOT):
    idx = rng.integers(0, n_val, size=n_val)
    yb = y_true_val[idx]
    sb = y_score_val[idx]
    mb = metrics_at_cutoff(yb, sb, cutoff)
    acc_b[b]  = mb["accuracy"]
    sens_b[b] = mb["sensitivity"]
    spec_b[b] = mb["specificity"]
    ppv_b[b]  = mb["ppv"]
    npv_b[b]  = mb["npv"]

acc_ci  = ci_from_boot(acc_b)
sens_ci = ci_from_boot(sens_b)
spec_ci = ci_from_boot(spec_b)
ppv_ci  = ci_from_boot(ppv_b)
npv_ci  = ci_from_boot(npv_b)

# --------------------
# Build summary table and save
# --------------------
rows = []

rows.append({
    "Dataset": "Training",
    "Cutoff_fraction": cutoff,
    "Cutoff_percent": cutoff * 100.0,
    "AUC": auc_tr,
    "TN": m_tr["tn"], "FP": m_tr["fp"], "FN": m_tr["fn"], "TP": m_tr["tp"],
    "Accuracy": m_tr["accuracy"],
    "Sensitivity": m_tr["sensitivity"],
    "Specificity": m_tr["specificity"],
    "PPV": m_tr["ppv"],
    "NPV": m_tr["npv"],
    "Accuracy_CI_low": np.nan, "Accuracy_CI_high": np.nan,
    "Sensitivity_CI_low": np.nan, "Sensitivity_CI_high": np.nan,
    "Specificity_CI_low": np.nan, "Specificity_CI_high": np.nan,
    "PPV_CI_low": np.nan, "PPV_CI_high": np.nan,
    "NPV_CI_low": np.nan, "NPV_CI_high": np.nan,
    "Cutoff_CI_low": cut_ci[0],
    "Cutoff_CI_high": cut_ci[1]
})

rows.append({
    "Dataset": "Validation",
    "Cutoff_fraction": cutoff,
    "Cutoff_percent": cutoff * 100.0,
    "AUC": auc_v,
    "TN": m_val["tn"], "FP": m_val["fp"], "FN": m_val["fn"], "TP": m_val["tp"],
    "Accuracy": m_val["accuracy"],
    "Sensitivity": m_val["sensitivity"],
    "Specificity": m_val["specificity"],
    "PPV": m_val["ppv"],
    "NPV": m_val["npv"],
    "Accuracy_CI_low": acc_ci[0], "Accuracy_CI_high": acc_ci[1],
    "Sensitivity_CI_low": sens_ci[0], "Sensitivity_CI_high": sens_ci[1],
    "Specificity_CI_low": spec_ci[0], "Specificity_CI_high": spec_ci[1],
    "PPV_CI_low": ppv_ci[0], "PPV_CI_high": ppv_ci[1],
    "NPV_CI_low": npv_ci[0], "NPV_CI_high": npv_ci[1],
    "Cutoff_CI_low": cut_ci[0],
    "Cutoff_CI_high": cut_ci[1]
})

summary_df = pd.DataFrame(rows)
summary_df.to_csv(OUT_METRICS, index=False)

# --------------------
# Export predictions with assigned classes (auditing)
# --------------------
train_out = train.copy()
train_out["Predicted_class"] = (train_out["Predictedlabel_fraction"] >= cutoff).astype(int)
train_out.to_csv(OUT_TRAIN_PRED, index=False)

val_out = val.copy()
val_out["Predicted_class"] = (val_out["Predictedlabel_fraction"] >= cutoff).astype(int)
val_out.to_csv(OUT_VAL_PRED, index=False)

print("\nSaved:")
print(f" - Metrics summary CSV     : {OUT_METRICS}")
print(f" - Training predictions CSV: {OUT_TRAIN_PRED}")
print(f" - Validation predictions  : {OUT_VAL_PRED}")
print(f" - Training ROC PNG        : {OUT_TRAIN_ROC}")
print(f" - Validation ROC PNG      : {OUT_VAL_ROC}")