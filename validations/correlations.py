# Created by alexandra at 13/06/2025
import pandas as pd
from scipy import stats

# 1. Load the data
df = pd.read_csv('data/classification_data.csv', delim_whitespace=True)

# 2. Continuous–continuous correlations between ALFAssay and ichorTF
pearson_ichor = df['ALFAssay'].corr(df['ichorTF'], method='pearson')
spearman_ichor = df['ALFAssay'].corr(df['ichorTF'], method='spearman')

# 3. Continuous–binary correlations between ALFAssay and oncoDNADetection
#    (Pearson with binary is a point-biserial; we also compute pointbiserialr directly)
pearson_onco = df['ALFAssay'].corr(df['oncoDNADetection'], method='pearson')
spearman_onco = df['ALFAssay'].corr(df['oncoDNADetection'], method='spearman')
pb_r, pb_p = stats.pointbiserialr(df['oncoDNADetection'], df['ALFAssay'])

# 4. Print the results
print(f"Correlation ALFAssay vs ichorTF:")
print(f"  Pearson   = {pearson_ichor:.3f}")
print(f"  Spearman  = {spearman_ichor:.3f}\n")

print(f"Correlation ALFAssay vs oncoDNADetection:")
print(f"  Pearson (point-biserial) = {pearson_onco:.3f}")
print(f"  Spearman                  = {spearman_onco:.3f}")
print(f"  pointbiserialr r, p-value = ({pb_r:.3f}, {pb_p:.3e})")




import pandas as pd
import numpy as np
from scipy import stats
from sklearn.metrics import roc_auc_score

# 1. Load the data
df = pd.read_csv('data/classification_data.csv', delim_whitespace=True)

# 2. Compute point-biserial correlation and p-value
r_pb, p_pb = stats.pointbiserialr(df['oncoDNADetection'], df['ALFAssay'])
print(f"Point-biserial r = {r_pb:.3f}, p-value = {p_pb:.3e}")

# 3. Interpret effect size using Cohen's benchmarks
r_abs = abs(r_pb)
if r_abs < 0.1:
    magnitude = "negligible"
elif r_abs < 0.3:
    magnitude = "small"
elif r_abs < 0.5:
    magnitude = "medium"
else:
    magnitude = "large"
print(f"Effect size: {magnitude} correlation")

# 4. Bootstrap CI for point-biserial r
n_boot = 1000
rng = np.random.RandomState(42)
boot_rs = []
for _ in range(n_boot):
    idxs = rng.choice(len(df), len(df), replace=True)
    r, _ = stats.pointbiserialr(df['oncoDNADetection'].values[idxs],
                                df['ALFAssay'].values[idxs])
    boot_rs.append(r)
ci_lower, ci_upper = np.percentile(boot_rs, [2.5, 97.5])
print(f"95% CI for r: [{ci_lower:.3f}, {ci_upper:.3f}]")

# 5. Assess predictive performance via AUC
auc_val = roc_auc_score(df['oncoDNADetection'], df['ALFAssay'])
print(f"ROC AUC for ALFAssay alone = {auc_val:.3f}")

# 6. Summary
print("\nSummary:")
print(f"- Point-biserial r = {r_pb:.3f} ({magnitude}, p={p_pb:.3e})")
print(f"- 95% CI for r = [{ci_lower:.3f}, {ci_upper:.3f}]")
print(f"- AUC = {auc_val:.3f}")

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats.mstats import winsorize

# 1. Load data
df = pd.read_csv('data/classification_data.csv', delim_whitespace=True)

# 2. Plot histogram of ichorTF
plt.figure()
plt.hist(df['ichorTF'].dropna(), bins=30)
plt.xlabel('ichorTF')
plt.ylabel('Frequency')
plt.title('Histogram of ichorTF')
plt.savefig("plots/histogram_of_ichorTF.png")
# plt.show()

# 3. Plot histogram of ALFAssay
plt.figure()
plt.hist(df['ALFAssay'].dropna(), bins=30)
plt.xlabel('ALFAssay')
plt.ylabel('Frequency')
plt.title('Histogram of ALFAssay')
plt.savefig("plots/histogram_of_ALFAssay.png")
# plt.show()

# 4. Boxplot of ichorTF
plt.figure()
plt.boxplot(df['ichorTF'].dropna(), vert=False)
plt.xlabel('ichorTF')
plt.title('Boxplot of ichorTF')
plt.savefig("plots/boxplot_of_ichorTF.png")
# plt.show()

# 5. Boxplot of ALFAssay
plt.figure()
plt.boxplot(df['ALFAssay'].dropna(), vert=False)
plt.xlabel('ALFAssay')
plt.title('Boxplot of ALFAssay')
plt.savefig("plots/boxplot_of_ALFAssay.png")
# plt.show()

# 6. Log-transform ALFAssay and plot histogram
df['log_ALFAssay'] = np.log1p(df['ALFAssay'])
plt.figure()
plt.hist(df['log_ALFAssay'].dropna(), bins=30)
plt.xlabel('log(ALFAssay + 1)')
plt.ylabel('Frequency')
plt.title('Histogram of log-transformed ALFAssay')
plt.savefig("plots/histogram_logT_of_ALFAssay.png")
# plt.show()

# 7. Winsorize ALFAssay at 1st and 99th percentiles and plot boxplot
df['alf_winsor'] = winsorize(df['ALFAssay'], limits=[0.01, 0.01])
plt.figure()
plt.boxplot(df['alf_winsor'].dropna(), vert=False)
plt.xlabel('Winsorized ALFAssay')
plt.title('Boxplot of Winsorized ALFAssay (1st–99th percentiles)')
plt.savefig("plots/boxplot_Winsorized_of_ALFAssay.png")
# plt.show()

