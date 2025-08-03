# Created by alexandra at 30/07/2025

import pandas as pd
import numpy as np

# Load input data
high_tf_samples = pd.read_csv('data/data_with_high_tf.csv', sep="\t")
healthy_samples = pd.read_csv('data/healthy_data.csv', sep="\t")
features = pd.read_csv('data/features.csv', index_col=0, sep="\t")
metadata = pd.read_csv('/Users/alexandra/PhD/FragmentationPatterns/Data/MetaData/AllMetaData_.csv', sep="\t")

# Calculate ctDNA concentration
metadata = metadata[['PatientId', 'ichorTF', 'VAF']].dropna()
metadata['VAF'] = metadata['VAF'] / 100
metadata['ctDNA_frac'] = metadata[['ichorTF', 'VAF']].mean(axis=1)

# Filter metadata to only high TF samples
high_tf_samples = high_tf_samples.merge(metadata[['PatientId', 'ctDNA_frac']], on='PatientId')

# Select corresponding features
healthy_features = features[healthy_samples['PatientId'].values]
tumor_features = features[high_tf_samples['PatientId'].values]

# Desired spike-in levels
spike_levels = [0.001, 0.005, 0.01, 0.05, 0.1]


# Create synthetic mixtures
def create_spike_in_sample(healthy, tumor, spike_fraction):
    return healthy * (1 - spike_fraction) + tumor * spike_fraction


synthetic_samples = {}
# Prepare list to collect synthetic metadata entries
synthetic_meta = []

for spike in spike_levels:
    synthetic_data = []
    synthetic_ids = []

    for healthy_id in healthy_features.columns:
        # Randomly pick one tumor sample
        tumor_id = np.random.choice(tumor_features.columns)

        healthy_profile = healthy_features[healthy_id]
        tumor_profile = tumor_features[tumor_id]

        synthetic_profile = create_spike_in_sample(healthy_profile, tumor_profile, spike)

        synthetic_id = f"synthetic_{healthy_id}_tumor_{tumor_id}_frac_{int(spike * 1000) / 10}%"
        synthetic_ids.append(synthetic_id)
        synthetic_data.append(synthetic_profile)

        # Record metadata for this synthetic sample
        synthetic_meta.append({
            'PatientId': synthetic_id,
            'ichorTF': spike,
            'VAF': spike,
            'ctDNADetected': 1,
            'Label': 1,
            'Process': 'Validation'
        })

    # Create dataframe and save
    synthetic_df = pd.DataFrame(synthetic_data).T
    synthetic_df.columns = synthetic_ids

    synthetic_samples[spike] = synthetic_df

    # Save each synthetic dataset
    synthetic_df.to_csv(f'data/synthetic_ctDNA_frac_{int(spike * 1000) / 10}pct.csv', sep="\t")

# Combine all into a single file (optional)
combined_synthetic = pd.concat(synthetic_samples.values(), axis=1)
combined_synthetic.to_csv('data/synthetic_ctDNA_combined_all_fractions.csv', sep="\t", index=False)

meta_df = pd.DataFrame(synthetic_meta)
meta_df.to_csv('data/synthetic_AllMetaData.csv', index=False,  sep="\t")

print("Synthetic datasets created and saved successfully.")

