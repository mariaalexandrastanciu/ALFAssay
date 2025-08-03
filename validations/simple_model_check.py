# Created by alexandra at 22/06/2025
import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import KFold, cross_val_predict

def create_data():
    features = pd.read_csv("/Users/alexandra/PhD/FragmentationPatterns/Data/NN_5MB/NN_5MB_pivot_with_ratios_.csv", sep="\t")
    features = features.fillna(0)

    labels_cols = ["PatientId","FragleTF", "Label",
    "ichorTF", "timepoint", "VAF", "OncoFollowDNA",
    "survivalStatus", "survivalTime", "VAFg0p001", "study",
    "MedianStatus", "batch", "Process", "status",
    "BCFS", "ichorCNA", "ctDNA", "Predictedlabel",
    "Truelabel", "ALFAssay", "Fragle", "OncoDNAVAF"]

    metadata = pd.read_csv("~/PhD/FragmentationPatterns/Data/MetaData/MetaDataAndPredictions.csv", sep="\t")
    metadata = metadata[metadata["PatientId"].notna()]
    metadata["ichorTF"] = metadata["ichorTF"].apply(float)
    metadata["ichorTF"] = np.where(metadata["ichorTF"] < 0.03, 0, metadata["ichorTF"] * 100)
    # metadata["Label"] = metadata["Label"].apply(float)
    metadata["OncoDNAVAF"] = metadata["OncoDNAVAF"].apply(float)
    metadata["Predictedlabel"] = metadata["Predictedlabel"].apply(float)
    metadata["FragleTF"] = metadata["FragleTF"].apply(float)
    patients = list(metadata["PatientId"].values)

    data_feat_cols = ["bin"] + ["Feature"] + patients
    data_feat = features[data_feat_cols]
    data_feat = data_feat[data_feat["Feature"]=="ratioShort"][patients].apply(pd.to_numeric)
    cols_feat = ["bin" + str(bin) for bin in list(range(1, 205, 1))]
    df_short = data_feat.T.reset_index()
    df_short.columns = ["PatientId"] + cols_feat
    df_short = df_short.merge(metadata, on="PatientId", how="left")
    df_short.to_csv("data_sorted_short_meta.csv", index=False, sep="\t")


def regression_model():

    # 1) Load your full data (including ichorTF + Process)
    df_train = pd.read_csv('../DataSummary/data_sorted_short_Training.csv', sep='\t', skiprows=[1,2])
    df_validation = pd.read_csv('../DataSummary/data_sorted_short_Validation.csv', sep='\t', skiprows=[1, 2])

    # 2) Compute median of bins 1..204 per sample
    bin_cols = [f"bin{i}" for i in range(1, 205)]
    df_train['median_bin'] = df_train[bin_cols].median(axis=1)
    df_validation['median_bin'] = df_validation[bin_cols].median(axis=1)

    # 3) Split into Training vs. Validation
    train = df_train.copy()
    val   = df_validation.copy()

    X_train, y_train = train[['median_bin']], train['ichorTF']
    X_val,   y_val   = val[['median_bin']],   val['ichorTF']

    # 4) Cross-validated predictions on Training set
    model = LinearRegression()
    kf = KFold(n_splits=5, shuffle=True, random_state=42)
    train_preds = cross_val_predict(model, X_train, y_train, cv=kf)

    # 5) Fit on all Training, predict Validation
    model.fit(X_train, y_train)
    val_preds = model.predict(X_val)

    # 6) Build results table
    train['prediction'] = train_preds
    train['Process'] = "Training"
    val['prediction']   = val_preds
    val['Process'] = "Validation"

    results = pd.concat([
        train[['PatientId', 'Process',  'ichorTF', 'prediction']],
        val  [['PatientId', 'Process', 'ichorTF', 'prediction']]
    ], ignore_index=True)

    print(results)
    results.to_csv("regression_results_short_t35.csv", index=False, sep="\t")

#
regression_model()
# create_data()