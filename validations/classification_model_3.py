import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from scipy.stats.mstats import winsorize
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler, PolynomialFeatures
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.model_selection import StratifiedKFold, GridSearchCV, cross_val_predict
from sklearn.metrics import roc_auc_score, roc_curve, auc

# 1. Load data
df = pd.read_csv('data/classification_data.csv', delim_whitespace=True)

# 2. Preprocess ALFAssay: winsorize at 1–99%, then log-transform
df['ALFA_winsor'] = winsorize(df['ALFAssay'], limits=[0.01, 0.01])
df['ALFA_wl'] = np.log1p(df['ALFA_winsor'])

# 3. Ground truth
y = df['oncoDNADetection'].astype(int).values

# 4. Scatter plot ichorTF vs ALFA_wl colored by label
plt.figure(figsize=(6,6))
plt.scatter(df['ichorTF'], df['ALFA_wl'], c=y, cmap='bwr', alpha=0.7)
plt.xlabel('ichorTF')
plt.ylabel('log1p(winsorized ALFAssay)')
plt.title('ichorTF vs ALFAssay')
plt.colorbar(label='oncoDNADetection')
plt.show()

# 5. Pearson correlation ichorTF vs ALFA_wl
r, p = pearsonr(df['ichorTF'], df['ALFA_wl'])
print(f"Correlation ichorTF vs ALFA_wl: r = {r:.3f}, p = {p:.3e}")

# 6. CV and grid setup
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
param_grid = {'svm__C': [0.1,1,10,100], 'svm__gamma': [0.01,0.1,1,10]}

# 7. Define pipelines
# Single feature
pipe_single = Pipeline([('scaler', StandardScaler()),
                        ('svm', SVC(kernel='rbf', probability=True))])
# Two features
pipe_two = Pipeline([('scaler', StandardScaler()),
                     ('svm', SVC(kernel='rbf', probability=True))])
# PCA + SVM
pipe_pca = Pipeline([('scale', StandardScaler()),
                     ('pca', PCA(n_components=2)),
                     ('svm', SVC(kernel='rbf', probability=True))])
# Poly features + SVM
pipe_poly = Pipeline([('scale', StandardScaler()),
                      ('poly', PolynomialFeatures(degree=2, include_bias=False)),
                      ('svm', SVC(kernel='rbf', probability=True))])
# L1 Logistic
pipe_l1 = Pipeline([('scale', StandardScaler()),
                    ('lr', LogisticRegression(penalty='l1', solver='saga', max_iter=5000))])

# 8. Feature matrices
X_single = df[['ichorTF']].values
X_two    = df[['ichorTF','ALFA_wl']].values
X_pca    = X_two.copy()  # PCA on these two
X_poly   = X_two.copy()
X_l1     = X_two.copy()

# 9. GridSearchCV for SVM pipelines
for name, pipe, X in [('Single', pipe_single, X_single),
                      ('Two',    pipe_two,    X_two),
                      ('PCA',    pipe_pca,    X_pca),
                      ('Poly',   pipe_poly,   X_poly)]:
    gs = GridSearchCV(pipe, param_grid, cv=cv, scoring='roc_auc')
    gs.fit(X, y)
    print(f"{name} SVM best params: {gs.best_params_}")
    # replace pipeline with best estimator
    if name=='Single': pipe_single = gs.best_estimator_
    if name=='Two':    pipe_two    = gs.best_estimator_
    if name=='PCA':    pipe_pca    = gs.best_estimator_
    if name=='Poly':   pipe_poly   = gs.best_estimator_

# 10. Cross-validated probabilities
probas = {}
for name, pipe, X in [('Single', pipe_single, X_single),
                      ('Two',    pipe_two,    X_two),
                      ('PCA',    pipe_pca,    X_pca),
                      ('Poly',   pipe_poly,   X_poly),
                      ('L1',     pipe_l1,     X_l1)]:
    probas[name] = cross_val_predict(pipe, X, y, cv=cv, method='predict_proba')[:,1]

# 11. Compute and display AUCs
aucs = {name: roc_auc_score(y, probas[name]) for name in probas}
for name, val in aucs.items():
    print(f"{name} model ROC AUC = {val:.3f}")

# 12. Paired bootstrap test: Two vs Single
def paired_bootstrap_auc(y, p1, p2, n_boot=1000):
    rng = np.random.RandomState(0)
    diffs = []
    for _ in range(n_boot):
        idx = rng.randint(0, len(y), len(y))
        if len(np.unique(y[idx]))<2:
            continue
        a1 = roc_auc_score(y[idx], p1[idx])
        a2 = roc_auc_score(y[idx], p2[idx])
        diffs.append(a2 - a1)
    lower, upper = np.percentile(diffs, [2.5,97.5])
    return lower, upper

low, up = paired_bootstrap_auc(y, probas['Single'], probas['Two'])
print(f"ΔAUC (Two–Single) 95% CI = [{low:.3f}, {up:.3f}]")
