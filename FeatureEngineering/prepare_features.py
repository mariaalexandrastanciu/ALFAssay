# Created by alexandra at 03/04/2025
import pandas as pd
from statsmodels.nonparametric.smoothers_lowess import lowess

from . import CONSTANTS as c
import glob
from . import utils as u
import numpy as np


def get_feature_columns():
    feature_columns = ["noReads", "noReadsGCParagonFrag", "ultraShortFrag", "ultraShortFragGCPara",
                       "shortFrag", "shortFragGCPara", "longFrag", "longFragGCPara", "ultraLongFrag",
                       "ultraLongFragGCPara", "shortOverLong", "shortOverLongGCPara"]

    return feature_columns


def concat_samples(input_files, output, meta_data, first_file, write_mode):
    i = 0
    for input_file in input_files:
        sample_name = u.get_patient_from_file_name(input_file, "_nn_data.bed")
        if sample_name in meta_data["PatientId"].values:
            if i > 0:
                first_file = False
                write_mode = "a"
            print("File processed: " + input_file)
            meta_data_per_patient = meta_data[meta_data["PatientId"] == sample_name]

            data = pd.read_csv(input_file, sep="\t")
            data["Chromosome"] = pd.Categorical(data['Chromosome'], c.chromosomes)
            data_sorted = data.sort_values(["Chromosome", "Start"])

            data_sorted.update(data_sorted[get_feature_columns()].fillna(0))
            data_sorted["PatientId"] = meta_data_per_patient["PatientId"].values[0]

            all_data = pd.merge(data_sorted, meta_data_per_patient, on="PatientId")
            all_data.to_csv(output, sep="\t", mode=write_mode, index=False, header=first_file)

            i = i + 1


def pivot_data(input_file, output_file, fields):
    data = pd.read_csv(input_file, sep="\t")

    patients = data["PatientId"].unique()
    pivot_data = pd.DataFrame(data["bin"].unique(), columns=["bin"])
    pivot_data.reset_index(drop=True, inplace=True)

    for i, field in enumerate(fields):
        pivot_data_per_field = pd.DataFrame(data["bin"].unique(), columns=["bin"])
        pivot_data.reset_index(drop=True, inplace=True)
        pivot_data_per_field["Feature"] = field
        for patient in patients:
            patient_data_reads = data[data["PatientId"] == patient][[field]]
            patient_data_reads = patient_data_reads.rename(columns={field: patient})
            patient_data_reads.reset_index(drop=True, inplace=True)
            pivot_data_per_field = pd.concat([pivot_data_per_field, patient_data_reads], axis=1)
        if i == 0:
            pivot_data = pivot_data_per_field
        else:
            pivot_data = pd.concat([pivot_data, pivot_data_per_field], axis=0)

    pivot_data.to_csv(output_file, index=False, sep="\t")


def loess_depth_residuals(bins_df: pd.DataFrame,
                          depth_df: pd.DataFrame,
                          frac: float = 0.3,
                          min_pts: int = 5) -> pd.DataFrame:
    """
    Remove non-linear depth-dependence from each bin via LOESS residuals,
    skipping bins with insufficient variability or points.

    Parameters
    ----------
    bins_df : pd.DataFrame
        Index = bin names,
        Columns = PatientId,
        Values = bin ratio for that patient.
    depth_df : pd.DataFrame
        Two columns: 'PatientId' and 'depth'.
    frac : float
        LOESS span parameter.
    min_pts : int
        Minimum number of non-NA points required to fit; otherwise skip.

    Returns
    -------
    resid_df : pd.DataFrame
        Same shape as bins_df, with (orig − loess_fit) or zero if skipped.
    """
    # 1) intersect patients
    common = bins_df.columns.intersection(depth_df['PatientId'])
    bins_sub = bins_df[common]
    depth_sub = depth_df.set_index('PatientId').loc[common, 'depth']

    # 2) prepare output
    resid_df = pd.DataFrame(0.0, index=bins_sub.index, columns=bins_sub.columns)

    x = depth_sub.values
    for bin_name in bins_sub.index:
        y = bins_sub.loc[bin_name].values

        # Identify valid points
        mask = ~np.isnan(y) & ~np.isnan(x)
        if mask.sum() < min_pts or np.nanstd(y[mask]) == 0:
            # Too few points or no variance: skip (residuals stay zero)
            continue

        try:
            # Fit LOESS on the valid subset
            fitted = lowess(endog=y[mask], exog=x[mask], frac=frac,
                            return_sorted=False)
            # Place residuals back into full array
            resid = np.full_like(y, np.nan, dtype=float)
            resid[mask] = y[mask] - fitted
            resid_df.loc[bin_name] = resid
        except Exception:
            # On any fitting error, fallback to mean-centering
            resid_df.loc[bin_name] = y - np.nanmean(y)

    return resid_df


def add_ratio():
    # 1. Load the pivoted table
    data = pd.read_csv('FeatureEngineering/data/NN_5MB_pivot.csv', sep="\t")
    metadata = pd.read_csv("refdata/AllMetaData.csv", sep="\t")
    patients = list(metadata["PatientId"].values)
    training_patients = metadata[metadata["Process"] == "Training"]["PatientId"].values
    validation_patients = metadata[metadata["Process"] == "Validation"]["PatientId"].values
    depth_df = metadata[["PatientId", "depth"]]
    data = data.fillna(0)
    eps = 0.001

    data_feat_reads = data[data["Feature"] == "noReadsGCParagonFrag"][patients].apply(pd.to_numeric) + eps
    data_feat_short_reads = data[data["Feature"] == "shortFragGCPara"][patients].apply(pd.to_numeric)
    data_feat_ultra_long_reads = data[data["Feature"] == "ultraLongFragGCPara"][patients].apply(pd.to_numeric)

    ratio_long = data_feat_ultra_long_reads.values / data_feat_reads
    ratio_short = data_feat_short_reads.values / data_feat_reads

    ##normalize short
    ratio_short_training = ratio_short[training_patients]
    ratio_short_validation = ratio_short[validation_patients]
    ratio_short_training_corr = loess_depth_residuals(ratio_short_training, depth_df)
    ratio_short_validation_corr = loess_depth_residuals(ratio_short_validation, depth_df)
    ratio_short_norm = pd.concat([ratio_short_training_corr, ratio_short_validation_corr], axis=1)

    ##normalize long
    ratio_long_training = ratio_long[training_patients]
    ratio_long_validation = ratio_long[validation_patients]
    ratio_long_training_corr = loess_depth_residuals(ratio_long_training, depth_df)
    ratio_long_validation_corr = loess_depth_residuals(ratio_long_validation, depth_df)
    ratio_long_norm = pd.concat([ratio_long_training_corr, ratio_long_validation_corr], axis=1)

    data2_short = pd.DataFrame(np.asarray(data[data["Feature"] == "noReadsGCParagonFrag"]["bin"].values).T,
                               columns=["bin"])
    data2_short["Feature"] = "ratioShort"
    # data2_short[patients] = ratio_short
    data2_short = pd.concat([data2_short, ratio_short_norm], axis=1)

    # data2_short_p = pd.DataFrame(  ratio_short_norm,  index=data2_short.index,   columns=[patients]   )
    # data2_short[patients] = data2_short_p

    data2_long = pd.DataFrame(np.asarray(data[data["Feature"] == "noReadsGCParagonFrag"]["bin"].values).T,
                              columns=["bin"])
    data2_long["Feature"] = "ratioLong"
    data2_long = pd.concat([data2_long, ratio_long_norm], axis=1)
    # data2_long_p = pd.DataFrame(  ratio_long_norm,  index=data2_long.index,   columns=[patients]   )
    # data2_long[patients] = data2_long_p

    full_data = pd.concat([data, data2_short, data2_long], axis=0)

    # 3. Compute the two new rows

    # 6. Save back out
    full_data.to_csv('FeatureEngineering/data/NN_5MB_features.csv', index=False, sep="\t")


def add_ratio_without_norm():
    # 1. Load the pivoted table
    data = pd.read_csv('FeatureEngineering/data/NN_5MB_pivot.csv', sep="\t")
    metadata = pd.read_csv("refdata/AllMetaData.csv", sep="\t")
    patients = list(metadata["PatientId"].values)
    data = data.fillna(0)
    eps = 0.001

    data_feat_reads = data[data["Feature"] == "noReadsGCParagonFrag"][patients].apply(pd.to_numeric) + eps
    data_feat_short_reads = data[data["Feature"] == "shortFragGCPara"][patients].apply(pd.to_numeric)
    data_feat_ultra_long_reads = data[data["Feature"] == "ultraLongFragGCPara"][patients].apply(pd.to_numeric)

    ratio_long = data_feat_ultra_long_reads.values / data_feat_reads
    ratio_short = data_feat_short_reads.values / data_feat_reads

    bin_values = np.asarray(data[data["Feature"] == "noReadsGCParagonFrag"]["bin"].values).T

    data2_short = pd.DataFrame(bin_values, columns=["bin"])
    data2_short["Feature"] = "ratioShort"
    data2_short = pd.concat([data2_short, ratio_short.reset_index(drop=True)], axis=1)

    data2_long = pd.DataFrame(bin_values, columns=["bin"])
    data2_long["Feature"] = "ratioLong"
    data2_long = pd.concat([data2_long, ratio_long.reset_index(drop=True)], axis=1)

    # Include noReadsGCParagonFrag rows from original data in the output
    data_reads_rows = data[data["Feature"] == "noReadsGCParagonFrag"].reset_index(drop=True)

    full_data = pd.concat([data, data2_short, data2_long, data_reads_rows], axis=0)

    # 6. Save back out
    full_data.to_csv('FeatureEngineering/data/NN_5MB_features_wn.csv', index=False, sep="\t")

def gather_features():
    first_file = True
    write_mode = "w"
    meta_data = pd.read_csv("refdata/AllMetaData.csv", sep="\t")

    input_files = sorted(glob.glob("FragmentsManipulations/data/output/*.bed"))
    output_folder = "FeatureEngineering/data/NN_5MB.csv"

    concat_samples(input_files, output_folder, meta_data, first_file, write_mode)

    fields = ["noReadsGCParagonFrag", "shortFragGCPara", "longFragGCPara", "shortOverLongGCPara", "ultraLongFragGCPara"]
    input_file = output_folder

    output_file = "FeatureEngineering/data/NN_5MB_pivot.csv"
    pivot_data(input_file, output_file, fields)
    add_ratio_without_norm()


gather_features()
