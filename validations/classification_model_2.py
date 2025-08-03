import pandas as pd
import numpy as np
from scipy.stats.mstats import winsorize
from sklearn.svm import SVC
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.model_selection import StratifiedKFold, GridSearchCV, cross_val_predict
from sklearn.metrics import roc_curve, auc, roc_auc_score
import matplotlib.pyplot as plt
import numpy as np
from sklearn.utils import check_array
from scipy import stats

def compute_midrank(x):
    """Computes midranks for tied values."""
    J = np.argsort(x)
    Z = x[J]
    N = len(x)
    T = np.zeros(N, dtype=float)
    i = 0
    while i < N:
        j = i
        while j < N and Z[j] == Z[i]:
            j += 1
        T[i:j] = 0.5*(i + j - 1)
        i = j
    T2 = np.empty(N, dtype=float)
    T2[J] = T + 1  # +1 for 1-based mid-ranks
    return T2

def fastDeLong(predictions_sorted_transposed, label_1_count):
    """
    Fast implementation of DeLong's algorithm for ROC AUC variance.
    Args:
        predictions_sorted_transposed: np.array shape (2, n)
        label_1_count: number of positive examples
    Returns:
        auc, auc_variance
    """
    m = label_1_count
    n = predictions_sorted_transposed.shape[1] - m
    tx = predictions_sorted_transposed[0, :]
    ty = predictions_sorted_transposed[1, :]
    vx = compute_midrank(tx)
    vy = compute_midrank(ty)
    auc = (vy[:m].sum() - m*(m+1)/2) / (m*n)
    v01 = (vy[:m] - vx[:m]) / n
    v10 = 1.0 - (vy[m:] - vx[m:]) / m
    sx = np.cov(v01, rowvar=False)
    sy = np.cov(v10, rowvar=False)
    s = sx/m + sy/n
    return auc, s

def delong_roc_test(y_true, y_pred1, y_pred2):
    """
    Computes p-value for two correlated ROC AUCs via DeLong test.
    """
    y_true = np.asarray(y_true)
    y_pred1 = np.asarray(y_pred1)
    y_pred2 = np.asarray(y_pred2)
    # sort by label then preds
    pos = y_pred1[y_true == 1]
    neg = y_pred1[y_true == 0]
    indices = np.argsort(np.concatenate([pos, neg]))
    label_1_count = len(pos)
    preds_sorted = np.vstack([np.concatenate([y_pred1[y_true == 1],
                                              y_pred1[y_true == 0]])[indices],
                              np.concatenate([y_pred2[y_true == 1],
                                              y_pred2[y_true == 0]])[indices]])
    auc1, var1 = fastDeLong(preds_sorted, label_1_count)
    auc2, var2 = fastDeLong(preds_sorted[::-1, :], label_1_count)  # swap predictions
    z = (auc1 - auc2) / np.sqrt(var1 + var2)
    p = stats.norm.sf(abs(z)) * 2
    return p


# 1. Load the data
df = pd.read_csv('data/classification_data.csv', delim_whitespace=True)

# 2. Winsorize ALFAssay at 1st and 99th percentiles
# df['ALFAssay_winsor'] = winsorize(df['ALFAssay'], limits=[0.01, 0.01])

# 3. Prepare feature matrices and target vector
X_ichor = df[['ichorTF']].values
X_both  = df[['ichorTF', 'ALFAssay']].values
X_alf   = df[['ALFAssay']].values
y       = df['oncoDNADetection'].astype(int).values

# 4. Stratified 5-fold CV
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

# 5. Base pipeline constructor
def make_pipeline():
    return Pipeline([
        ('scaler', StandardScaler()),
        ('svm',    SVC(kernel='rbf', probability=True, random_state=42))
    ])

# 6. Hyperparameter grid
param_grid = {
    'svm__C':     [0.1, 1, 10, 100],
    'svm__gamma': [0.01, 0.1, 1, 10]
}

# 7. Grid search for each model
grid1 = GridSearchCV(make_pipeline(), param_grid, cv=cv, scoring='roc_auc')
grid1.fit(X_ichor, y)
best1 = grid1.best_estimator_

grid2 = GridSearchCV(make_pipeline(), param_grid, cv=cv, scoring='roc_auc')
grid2.fit(X_both, y)
best2 = grid2.best_estimator_

grid3 = GridSearchCV(make_pipeline(), param_grid, cv=cv, scoring='roc_auc')
grid3.fit(X_alf, y)
best3 = grid3.best_estimator_

# 8. Out-of-fold probability predictions using tuned pipelines
proba1 = cross_val_predict(best1, X_ichor, y, cv=cv, method='predict_proba')[:, 1]
proba2 = cross_val_predict(best2, X_both,  y, cv=cv, method='predict_proba')[:, 1]
proba3 = cross_val_predict(best3, X_alf,   y, cv=cv, method='predict_proba')[:, 1]

# 8b. Flip probabilities if AUC < 0.5
for name, proba in [('ichorTF', proba1), ('ichorTF+ALFAssay', proba2), ('ALFAssay', proba3)]:
    if roc_auc_score(y, proba) < 0.5:
        if name == 'ichorTF':
            proba1 = 1 - proba1
        elif name == 'ichorTF+ALFAssay':
            proba2 = 1 - proba2
        else:
            proba3 = 1 - proba3

p_value = delong_roc_test(y, proba1, proba2)
print(f"DeLong test p-value between ichorTF and ichorTF+ALFAssay: {p_value:.3e}")
# 9. Compute ROC curves and AUCs
fpr1, tpr1, _ = roc_curve(y, proba1); auc1 = auc(fpr1, tpr1)
fpr2, tpr2, _ = roc_curve(y, proba2); auc2 = auc(fpr2, tpr2)
fpr3, tpr3, _ = roc_curve(y, proba3); auc3 = auc(fpr3, tpr3)

# 10. Bootstrap for 95% CI bands and AUC CIs
n_bootstraps = 1000
rng = np.random.RandomState(42)
fpr_grid = np.linspace(0, 1, 100)
tprs1, tprs2, tprs3 = [], [], []
auc1_boot, auc2_boot, auc3_boot = [], [], []

for _ in range(n_bootstraps):
    idxs = rng.randint(0, len(y), len(y))
    if len(np.unique(y[idxs])) < 2:
        continue
    fpr_bs1, tpr_bs1, _ = roc_curve(y[idxs], proba1[idxs])
    fpr_bs2, tpr_bs2, _ = roc_curve(y[idxs], proba2[idxs])
    fpr_bs3, tpr_bs3, _ = roc_curve(y[idxs], proba3[idxs])
    tprs1.append(np.interp(fpr_grid, fpr_bs1, tpr_bs1))
    tprs2.append(np.interp(fpr_grid, fpr_bs2, tpr_bs2))
    tprs3.append(np.interp(fpr_grid, fpr_bs3, tpr_bs3))
    auc1_boot.append(roc_auc_score(y[idxs], proba1[idxs]))
    auc2_boot.append(roc_auc_score(y[idxs], proba2[idxs]))
    auc3_boot.append(roc_auc_score(y[idxs], proba3[idxs]))

tprs1, tprs2, tprs3 = map(np.array, (tprs1, tprs2, tprs3))
tpr1_lower, tpr1_upper = np.percentile(tprs1, [2.5, 97.5], axis=0)
tpr2_lower, tpr2_upper = np.percentile(tprs2, [2.5, 97.5], axis=0)
tpr3_lower, tpr3_upper = np.percentile(tprs3, [2.5, 97.5], axis=0)
auc1_lower, auc1_upper = np.percentile(auc1_boot, [2.5, 97.5])
auc2_lower, auc2_upper = np.percentile(auc2_boot, [2.5, 97.5])
auc3_lower, auc3_upper = np.percentile(auc3_boot, [2.5, 97.5])

# 11. Plot ROC curves with CI bands and AUC [CI]
plt.figure(figsize=(6, 6))
plt.fill_between(fpr_grid, tpr1_lower, tpr1_upper, color='blue', alpha=0.2)
plt.plot(fpr1, tpr1, color='blue', lw=2,
         label=f'ichorTF only (AUC = {auc1:.3f} [{auc1_lower:.3f}-{auc1_upper:.3f}])')
plt.fill_between(fpr_grid, tpr2_lower, tpr2_upper, color='red', alpha=0.2)
plt.plot(fpr2, tpr2, color='red', lw=2,
         label=f'ichorTF+ALFAssay (AUC = {auc2:.3f} [{auc2_lower:.3f}-{auc2_upper:.3f}])')
plt.fill_between(fpr_grid, tpr3_lower, tpr3_upper, color='green', alpha=0.2)
plt.plot(fpr3, tpr3, color='green', lw=2,
         label=f'ALFAssay only (AUC = {auc3:.3f} [{auc3_lower:.3f}-{auc3_upper:.3f}])')
plt.plot([0, 1], [0, 1], color='grey', lw=1, linestyle='--')
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
# plt.title('Tuned SVM ROC Curves')
plt.legend(loc='lower right')
plt.axis('equal')
plt.tight_layout()
plt.savefig("plots/SVM_classifier_tuned.pdf")

# 12. Assemble per-patient predictions
results = pd.DataFrame({
    'PatientId':              df['PatientId'],
    'TrueLabel':              y,
    'Prob_ichorTF':           proba1,
    'Pred_ichorTF':           (proba1 >= 0.5).astype(int),
    'Prob_ichorTF_ALFAssay':  proba2,
    'Pred_ichorTF_ALFAssay':  (proba2 >= 0.5).astype(int),
    'Prob_ALFAssay':          proba3,
    'Pred_ALFAssay':          (proba3 >= 0.5).astype(int)
})

# 13. Statistical test for ΔAUC between ichorTF and ichorTF+ALFAssay
def paired_bootstrap_auc(y_true, p1, p2, n_boot=1000):
    rng = np.random.RandomState(0)
    diffs = []
    for _ in range(n_boot):
        idx = rng.randint(0, len(y_true), len(y_true))
        if len(np.unique(y_true[idx])) < 2:
            continue
        a1 = roc_auc_score(y_true[idx], p1[idx])
        a2 = roc_auc_score(y_true[idx], p2[idx])
        diffs.append(a2 - a1)
    return np.percentile(diffs, [2.5, 97.5]), np.mean(diffs)

(ci_low, ci_high), delta_mean = paired_bootstrap_auc(y, proba1, proba2)
print(f"ΔAUC (ichorTF+ALFAssay - ichorTF): mean = {delta_mean:.3f}, 95% CI = [{ci_low:.3f}, {ci_high:.3f}]")

print(results.head())

results.to_csv("results/classifier_predictions.csv", index=False, sep="\t")



# Example usage:

