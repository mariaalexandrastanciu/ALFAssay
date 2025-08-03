# Created by alexandra at 13/06/2025
import pandas as pd
import numpy as np
from scipy.stats.mstats import winsorize
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.model_selection import StratifiedKFold, GridSearchCV, cross_val_predict
from sklearn.metrics import roc_curve, auc, roc_auc_score
import matplotlib.pyplot as plt

# 1. Load the data
df = pd.read_csv('data/classification_data.csv', delim_whitespace=True)

# 2. Winsorize ALFAssay at 1st and 99th percentiles
df['ALFAssay_winsor'] = winsorize(df['ALFAssay'], limits=[0.01, 0.01])

# 3. Prepare feature matrices and target vector
X_ichor = df[['ichorTF']].values
X_both  = df[['ichorTF', 'ALFAssay_winsor']].values
X_alf   = df[['ALFAssay_winsor']].values
y       = df['oncoDNADetection'].astype(int).values

# 4. Stratified 5-fold CV
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

# 5. Base pipeline constructor for L1-penalized logistic regression
def make_l1_pipeline():
    return Pipeline([
        ('scaler', StandardScaler()),
        ('lr',     LogisticRegression(penalty='l1', solver='saga', max_iter=5000))
    ])

# 6. Hyperparameter grid for logistic regression (inverse regularization strength)
param_grid = {
    'lr__C': [0.01, 0.1, 1, 10, 100]
}

# 7. Grid search for each model
grid1 = GridSearchCV(make_l1_pipeline(), param_grid, cv=cv, scoring='roc_auc')
grid1.fit(X_ichor, y)
best1 = grid1.best_estimator_
print("Best params (ichorTF only):", grid1.best_params_)

grid2 = GridSearchCV(make_l1_pipeline(), param_grid, cv=cv, scoring='roc_auc')
grid2.fit(X_both, y)
best2 = grid2.best_estimator_
print("Best params (ichorTF + ALFAssay):", grid2.best_params_)

grid3 = GridSearchCV(make_l1_pipeline(), param_grid, cv=cv, scoring='roc_auc')
grid3.fit(X_alf, y)
best3 = grid3.best_estimator_
print("Best params (ALFAssay only):", grid3.best_params_)

# 8. Out-of-fold probability predictions using tuned pipelines
proba1 = cross_val_predict(best1, X_ichor, y, cv=cv, method='predict_proba')[:, 1]
proba2 = cross_val_predict(best2, X_both,  y, cv=cv, method='predict_proba')[:, 1]
proba3 = cross_val_predict(best3, X_alf,   y, cv=cv, method='predict_proba')[:, 1]

# 9. Flip probabilities if AUC < 0.5
for name, proba in [('ichorTF', proba1), ('ichorTF+ALFAssay', proba2), ('ALFAssay', proba3)]:
    if roc_auc_score(y, proba) < 0.5:
        if name == 'ichorTF':
            proba1 = 1 - proba1
        elif name == 'ichorTF+ALFAssay':
            proba2 = 1 - proba2
        else:
            proba3 = 1 - proba3

# 10. Compute ROC curves and AUCs
fpr1, tpr1, _ = roc_curve(y, proba1); auc1 = auc(fpr1, tpr1)
fpr2, tpr2, _ = roc_curve(y, proba2); auc2 = auc(fpr2, tpr2)
fpr3, tpr3, _ = roc_curve(y, proba3); auc3 = auc(fpr3, tpr3)

# 11. Bootstrap for 95% CI bands and AUC CIs
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

# 12. Plot ROC curves with CI bands and AUC [CI]
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
plt.title('L1-penalized Logistic Regression ROC Curves')
plt.legend(loc='lower right')
plt.axis('equal')
plt.tight_layout()
plt.savefig("plots/L1_logistic_classifier.png")

# 13. Assemble per-patient predictions
results = pd.DataFrame({
    'PatientId':            df['PatientId'],
    'TrueLabel':            y,
    'Prob_ichorTF':         proba1,
    'Pred_ichorTF':         (proba1 >= 0.5).astype(int),
    'Prob_ichorTF_ALFAssay':proba2,
    'Pred_ichorTF_ALFAssay':(proba2 >= 0.5).astype(int),
    'Prob_ALFAssay':        proba3,
    'Pred_ALFAssay':        (proba3 >= 0.5).astype(int)
})

print(results.head())
