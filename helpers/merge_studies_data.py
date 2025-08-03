# Created by alexandra at 03/04/2025
import pandas as pd
import CONSTANTS as c
import glob
import utils as u
from scipy.ndimage import gaussian_filter1d
import numpy as np

def get_feature_columns():
    feature_columns = ["noReads", "noReadsGCParagonFrag", "ultraShortFrag", "ultraShortFragGCPara",
                       "shortFrag", "shortFragGCPara", "longFrag", "longFragGCPara", "ultraLongFrag",
                       "ultraLongFragGCPara", "shortOverLong", "shortOverLongGCPara"]

    return feature_columns

def concat_samples(input_files, output, meta_data, study, first_file, write_mode):
    i=0
    for input_file in input_files:
        sample_name = u.get_patient_from_file_name(input_file, "_nn_data.bed")  # os.path.basename(input_file).replace("_fragment_size_summary_window.csv", "")
        # if study == "healthy_sWGS":
        #     sample_name = sample_name.split("-xx", 1)[0]
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
            all_data["study"] = study
            all_data.to_csv(output, sep="\t", mode=write_mode, index=False, header=first_file)

            i = i + 1


def pivot_data(input_file, output_file, fields):
    data = pd.read_csv(input_file, sep ="\t")

    patients = data["PatientId"].unique()
    pivot_data = pd.DataFrame(data["bin"].unique(), columns=["bin"])
    pivot_data.reset_index(drop=True, inplace=True)

    for i, field in enumerate(fields):
        pivot_data_per_field = pd.DataFrame(data["bin"].unique(), columns=["bin"])
        pivot_data.reset_index(drop=True, inplace=True)
        pivot_data_per_field["Feature"] = field
        for patient in patients:
            patient_data_reads = data[data["PatientId"]==patient][[field]]
            patient_data_reads = patient_data_reads.rename(columns={field: patient})
            patient_data_reads.reset_index(drop=True, inplace=True)
            pivot_data_per_field = pd.concat([pivot_data_per_field, patient_data_reads], axis=1)
        if i == 0:
            pivot_data = pivot_data_per_field
        else:
            pivot_data = pd.concat([pivot_data, pivot_data_per_field], axis=0)

    pivot_data.to_csv(output_file, index=False, sep="\t")

def run_input_file():
    studies = ["healthy_sWGS", "NeoRheaStudy", "PearlStudy", "SynergyStudy"]
    root = "/Users/alexandra/PhD/"
    first_file = True
    write_mode = "w"
    analysis_input = "NN_5MB" #"WindowFragmentSizes_30_700_mq40_wz100k" #"WindowFragmentSizes_hg38_30_700_mq60" #"WindowFragmentSizes_30_700_mq40_wz100k"
    analysis_output = "NN_5MB"
    # meta_data = pd.read_csv("/Users/alexandra/PhD/FragmentationPatterns/Data/MetaData/AllStudiesMetaDataWDELFI.csv", sep="\t")
    for i, study in enumerate(studies):
        if study == "healthy_sWGS":
            meta = "Healthy"
        else:
            meta = study
        if i > 0:
            first_file = False
            write_mode = "a"
        meta_data = pd.read_csv(root + "FragmentationPatterns/Data/MetaData/" + meta + "MetaData.csv", sep="\t")
        # meta_data = meta_data[meta_data["Label"] == "Healthy_Synergy"]

        input_files = sorted(glob.glob(root + study +
                                "/FragmentationPatterns/FragmentSizes/" + analysis_input + "/*.bed"))
        output_folder = root + "FragmentationPatterns/Data/" + analysis_output + "/NN_5MB.csv"

        concat_samples(input_files, output_folder, meta_data, study, first_file, write_mode)

    fields = ["noReadsGCParagonFrag", "shortFragGCPara", "longFragGCPara", "shortOverLongGCPara", "ultraLongFragGCPara"]
    root = "/Users/alexandra/PhD/FragmentationPatterns/Data/NN_5MB/"
    input_file = root + "NN_5MB.csv"

    output_file = root + "NN_5MB_pivot.csv"
    pivot_data(input_file, output_file, fields)

# run_input_file()

def smooth_data(pivot_file):
    # 1. Load your pivoted data (bins × samples)
    data = pd.read_csv(pivot_file, sep='\t')

    meta_cols = ['bin', 'Feature']
    sample_cols = [c for c in data.columns if c not in meta_cols]

    # 3. Extract the numeric matrix and index by bin
    df_bins = data[sample_cols].astype(float)
    df_bins.index = data['bin']

    # 4. Compute a tiny pseudocount ε = half the smallest non-zero entry
    min_nonzero = df_bins[df_bins > 0].min().min()
    epsilon = min_nonzero * 0.5

    # 5. Add pseudocount to all values
    df_pc = df_bins + epsilon

    # 6. Mask the original zeros (so only those get interpolated)
    df_masked = df_pc.mask(df_bins == 0, np.nan)

    # 7. Linearly interpolate down each sample column (across bins)
    df_interp = df_masked.interpolate(
        axis=0,
        method='linear',
        limit_direction='both'
    )

    # 8. Fill any remaining NaNs at the ends by forward/backward fill
    df_filled = df_interp.apply(
        lambda col: col.fillna(method='ffill').fillna(method='bfill'),
        axis=0
    )

    # 9. Smooth each sample’s profile across bins with a Gaussian filter (σ = 2)
    smoothed_arr = gaussian_filter1d(df_filled.values, sigma=2, axis=0)
    df_smoothed = pd.DataFrame(
        smoothed_arr,
        index=df_filled.index,
        columns=sample_cols
    )

    # 10. Reassemble into a full DataFrame with original column order
    df_corrected_full = pd.concat([
        data[['bin', 'Feature']].reset_index(drop=True),
        df_smoothed.reset_index(drop=True)
    ], axis=1)[data.columns]

    df_corrected_full.to_csv("/Users/alexandra/PhD/FragmentationPatterns/Data/NN_5MB/NN_5MB_pivot_v1.csv", index=False, sep="\t")


# pivot_file ="/Users/alexandra/PhD/FragmentationPatterns/Data/NN_5MB/NN_5MB_pivot.csv"
# smooth_data(pivot_file)


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
    data = pd.read_csv('/Users/alexandra/PhD/FragmentationPatterns/Data/NN_5MB/NN_5MB_pivot.csv', sep="\t")
    metadata = pd.read_csv("/Users/alexandra/PhD/FragmentationPatterns/Data/MetaData/AllMetaData_.csv", sep="\t")
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
    ratio_short_norm = pd.concat([ratio_short_training_corr,ratio_short_validation_corr],axis=1 )

    ##normalize long
    ratio_long_training = ratio_long[training_patients]
    ratio_long_validation = ratio_long[validation_patients]
    ratio_long_training_corr = loess_depth_residuals(ratio_long_training, depth_df)
    ratio_long_validation_corr = loess_depth_residuals(ratio_long_validation, depth_df)
    ratio_long_norm = pd.concat([ratio_long_training_corr, ratio_long_validation_corr], axis=1)

    data2_short = pd.DataFrame(np.asarray(data[data["Feature"] == "noReadsGCParagonFrag"]["bin"].values).T, columns=["bin"] )
    data2_short["Feature"] = "ratioShort"
    # data2_short[patients] = ratio_short
    data2_short = pd.concat([data2_short,ratio_short_norm],axis=1 )

    # data2_short_p = pd.DataFrame(  ratio_short_norm,  index=data2_short.index,   columns=[patients]   )
    # data2_short[patients] = data2_short_p

    data2_long = pd.DataFrame(np.asarray(data[data["Feature"] == "noReadsGCParagonFrag"]["bin"].values).T, columns=["bin"] )
    data2_long["Feature"] = "ratioLong"
    data2_long = pd.concat([data2_long, ratio_long_norm], axis=1)
    # data2_long_p = pd.DataFrame(  ratio_long_norm,  index=data2_long.index,   columns=[patients]   )
    # data2_long[patients] = data2_long_p


    full_data =   pd.concat([data, data2_short, data2_long], axis=0)

    # 3. Compute the two new rows

    # 6. Save back out
    full_data.to_csv('/Users/alexandra/PhD/FragmentationPatterns/Data/NN_5MB/NN_5MB_pivot_with_ratios_.csv', index=False, sep="\t")
    # print("Saved augmented table with 'ratioShort' and 'ratioLong' to NN_5MB_pivot_with_ratios.csv")

add_ratio()
