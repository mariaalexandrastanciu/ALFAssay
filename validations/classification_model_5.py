import pandas as pd
from sklearn.svm import SVC
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.model_selection import StratifiedKFold, GridSearchCV, cross_val_predict
from sklearn.metrics import roc_curve, auc, roc_auc_score
import matplotlib.pyplot as plt

# -- DeLong functions (unchanged) --
import numpy as np
import scipy.stats as st

def compute_midrank(x):
    J = np.argsort(x)
    Z = x[J]
    N = len(x)
    T = np.zeros(N, float)
    i = 0
    while i < N:
        j = i
        while j < N and Z[j] == Z[i]:
            j += 1
        T[i:j] = 0.5 * (i + j - 1)
        i = j
    T2 = np.empty(N, float)
    T2[J] = T + 1
    return T2

def fast_delong(preds, m):
    n = preds.shape[1] - m
    tx, ty = preds
    vx = compute_midrank(tx)
    vy = compute_midrank(ty)
    auc = (vy[:m].sum() - m*(m+1)/2) / (m*n)
    v01 = (vy[:m] - vx[:m]) / n
    v10 = 1 - (vy[m:] - vx[m:]) / m
    sx = np.cov(v01, rowvar=False)
    sy = np.cov(v10, rowvar=False)
    return auc, sx/m + sy/n

def delong_roc_test(y_true, p1, p2):
    y_true = np.asarray(y_true)
    pos = p1[y_true==1]; neg = p1[y_true==0]
    m = len(pos)
    preds = np.vstack([
        np.concatenate([pos, neg]),
        np.concatenate([p2[y_true==1], p2[y_true==0]])
    ])
    auc1, var1 = fast_delong(preds, m)
    auc2, var2 = fast_delong(preds[::-1], m)
    z = (auc1 - auc2) / np.sqrt(var1 + var2)
    p = st.norm.sf(abs(z))*2
    return p

# 1. Load data
df = pd.read_csv('data/classification_data_ichorTF.csv', delim_whitespace=True)

# 2. Ground truth: ichorTF (binary 0/1)
y = df['ichorTF'].astype(int).values

# 3. Features: oncoDNA detection and ALFAssay
X_onco = df[['oncoDNADetection']].astype(int).values
X_alf  = df[['ALFAssay']].values
X_both = df[['oncoDNADetection', 'ALFAssay']].values

# 4. CV and pipeline
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
def make_pipeline():
    return Pipeline([
        ('scaler', StandardScaler()),
        ('svm',    SVC(kernel='rbf', probability=True, random_state=42))
    ])

param_grid = {
    'svm__C':     [0.1, 1, 10, 100],
    'svm__gamma': [0.01, 0.1, 1, 10]
}

# 5. GridSearch for each feature set
grid_onco = GridSearchCV(make_pipeline(), param_grid, cv=cv, scoring='roc_auc')
grid_onco.fit(X_onco, y)
best_onco = grid_onco.best_estimator_

grid_alf = GridSearchCV(make_pipeline(), param_grid, cv=cv, scoring='roc_auc')
grid_alf.fit(X_alf, y)
best_alf = grid_alf.best_estimator_

grid_both = GridSearchCV(make_pipeline(), param_grid, cv=cv, scoring='roc_auc')
grid_both.fit(X_both, y)
best_both = grid_both.best_estimator_

# 6. Out-of-fold probabilities
p_onco = cross_val_predict(best_onco, X_onco, y, cv=cv, method='predict_proba')[:,1]
p_alf  = cross_val_predict(best_alf,  X_alf,  y, cv=cv, method='predict_proba')[:,1]
p_both = cross_val_predict(best_both, X_both, y, cv=cv, method='predict_proba')[:,1]

# 6b. Flip if inverted
for name, p in [('oncoDNA', p_onco), ('ALFA', p_alf), ('both', p_both)]:
    if roc_auc_score(y, p) < 0.5:
        if name=='oncoDNA': p_onco = 1-p_onco
        if name=='ALFA':    p_alf  = 1-p_alf
        if name=='both':    p_both = 1-p_both

# 7. ROC and AUC
fpr_o, tpr_o, _ = roc_curve(y, p_onco);  auc_o  = auc(fpr_o, tpr_o)
fpr_a, tpr_a, _ = roc_curve(y, p_alf);   auc_a  = auc(fpr_a, tpr_a)
fpr_b, tpr_b, _ = roc_curve(y, p_both);  auc_b  = auc(fpr_b, tpr_b)

# 8. DeLong p-value
p_val = delong_roc_test(y, p_onco, p_both)

# 9. Plot
plt.figure(figsize=(6,6))
plt.plot(fpr_o, tpr_o, 'b-', label=f'oncoDNA only (AUC={auc_o:.3f})')
plt.plot(fpr_a, tpr_a, 'g-', label=f'ALFA only (AUC={auc_a:.3f})')
plt.plot(fpr_b, tpr_b, 'r-', label=f'both (AUC={auc_b:.3f})')
plt.plot([0,1],[0,1],'k--')
plt.xlabel('False Positive Rate'); plt.ylabel('True Positive Rate')
plt.legend(loc='lower right')
plt.title(f'ROC (p={p_val:.3e} by DeLong)')
plt.axis('equal'); plt.tight_layout()
plt.savefig("plots/SVM_classifier_tuned_ichorTF.pdf")
# plt.show()

# 10. Predictions table
results = pd.DataFrame({
    'PatientId':             df['PatientId'],
    'True_ichorTF':          y,
    'Prob_oncoDNA':          p_onco,
    'Prob_ALFA':             p_alf,
    'Prob_both':             p_both
})
results.to_csv("classifier_predictions_ichorTF_as_truth.tsv", sep='\t', index=False)
