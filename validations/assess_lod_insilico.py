# Created by alexandra at 01/08/2025
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc, classification_report, confusion_matrix
import statsmodels.api as sm

# 1. Load model results and metadata
results = pd.read_csv('/Users/alexandra/PhD/PyCharmProjects/ALFAssayNN/results/results_validationInSilico_ichorTF_t35.csv', sep="\t")
meta = pd.read_csv('data/synthetic_AllMetaData.csv', sep="\t")
# Merge on PatientId
df = results.merge(meta, on='PatientId')

# 2. Adjust true and predicted fractions
# True label is in percent: convert to fraction
df['TrueFrac'] = df['Truelabel'] / 100.0
# Predicted label: set negatives to zero
df['PredScore'] = df['Predictedlabel'].clip(lower=0)

# 3. Define binary detection based on 3% true fraction threshold
df['TrueDetect'] = (df['TrueFrac'] >= 0.03).astype(int)
# Predicted detection at same threshold
df['PredDetect'] = (df['PredScore'] >= 0.03).astype(int)

# 4. ROC curve & AUC for detection
fpr, tpr, thresholds = roc_curve(df['TrueDetect'], df['PredScore'])
roc_auc = auc(fpr, tpr)
print(f'ROC AUC: {roc_auc:.3f}')

# Plot ROC curve
plt.figure()
plt.plot(fpr, tpr, lw=2, label=f'AUC = {roc_auc:.3f}')
plt.plot([0, 1], [0, 1], '--', lw=1)
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('ROC Curve for ctDNA Detection')
plt.legend(loc='lower right')
plt.show()

# 5. Classification report
print(classification_report(df['TrueDetect'], df['PredDetect'],
                            target_names=['Negative', 'Positive']))

# 6. Sensitivity/Specificity by true fraction
rates = []
for frac in sorted(df['TrueFrac'].unique()):
    sub = df[df['TrueFrac'] == frac]
    tn, fp, fn, tp = confusion_matrix(sub['TrueDetect'], sub['PredDetect']).ravel()
    sens = tp / (tp + fn) if (tp + fn) > 0 else np.nan
    spec = tn / (tn + fp) if (tn + fp) > 0 else np.nan
    rates.append({'TrueFrac': frac, 'Sensitivity': sens, 'Specificity': spec})
rates_df = pd.DataFrame(rates)
print(rates_df)

# 7. Estimate Limit of Detection (LOD) at 95%
# Logistic regression on binary detection vs true fraction
X = sm.add_constant(df['TrueFrac'])
model = sm.GLM(df['PredDetect'], X, family=sm.families.Binomial()).fit()
# Solve for P(detected)=0.95
target_p = 0.95
logit_val = np.log(target_p / (1 - target_p))
intercept, slope = model.params['const'], model.params['TrueFrac']
lod = (logit_val - intercept) / slope
print(f'Estimated LOD at 95% detection: {lod:.4f}')

# 8. Plot detection probability curve
frac_range = np.linspace(df['TrueFrac'].min(), df['TrueFrac'].max(), 100)
X_pred = sm.add_constant(frac_range)
pred_prob = model.predict(X_pred)

plt.figure()
plt.plot(frac_range, pred_prob, lw=2)
plt.axhline(0.95, color='red', linestyle='--', label='95% prob')
plt.axvline(lod, color='green', linestyle='--', label=f'LOD = {lod:.4f}')
plt.xlabel('True ctDNA Fraction')
plt.ylabel('Detection Probability')
plt.title('Detection Probability vs True ctDNA Fraction')
plt.legend()
plt.show()