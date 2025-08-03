# Created by alexandra at 15/11/2024
import numpy as np
import pandas as pd
import pyranges as pr
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.pyplot import gcf
from scipy.stats import kstest
from scipy.ndimage import gaussian_filter1d
import glob
import os.path
from NNModel import utils as u
from statsmodels.nonparametric.smoothers_lowess import lowess

import matplotlib.pyplot as plt
from scipy.stats import linregress

# 1. assume bins_df and depth_df are already loaded:
#    bins_df has columns ['PatientId', 'bin1', ..., 'bin204']
#    depth_df has two unnamed columns: first is filename+PatientId, second is depth

# Rename depth_df columns and clean PatientId


def corr_depth_loess(bins_df,metadata, process):
    # depth_df = pd.read_csv("/Users/alexandra/PhD/healthy_sWGS/mosdepth/total_depth_healthy.txt", sep=" ", header=None)
    # depth = (
    #     depth_df
    #         .copy()
    #         .rename(columns={depth_df.columns[0]: 'raw_id', depth_df.columns[1]: 'depth'})
    # )
    # depth['PatientId'] = depth['raw_id'].str.replace(
    #     r'\.mosdepth\.summary\.txt$', '', regex=True
    # )
    # depth = depth[['PatientId', 'depth']]
    #
    # # 2) Merge into one DataFrame
    merged = bins_df.merge(metadata, on='PatientId', how='inner')

    # 3) LOESS‑fit each bin vs. depth, compute residuals
    bin_cols = [f'bin{i}' for i in range(1, 205)]
    # frac = 0.3  # adjust span as needed

    # resid_df = pd.DataFrame(index=merged.index, columns=bin_cols, dtype=float)
    #
    # for col in bin_cols:
    #     y = merged[col].values
    #     x = merged['depth'].values
    #     # fit LOESS: returns smoothed y-values in original order
    #     fitted = lowess(endog=y, exog=x, frac=frac, return_sorted=False)
    #     resid_df[col] = y - fitted
    #
    # # 4) Build corrected DataFrame
    # corrected = pd.concat([
    #     merged[['PatientId', 'depth']].reset_index(drop=True),
    #     resid_df.reset_index(drop=True)
    # ], axis=1)

    # 5) Compute per-patient median of LOESS residuals
    merged['median_loess_resid'] = merged[bin_cols].median(axis=1)

    # 6) Check correlation and plot
    r_loess = merged['median_loess_resid'].corr(merged['depth'])
    print(f'LOESS‑residual median vs. depth Pearson r = {r_loess:.3f}')

    plt.figure(figsize=(6, 6))
    plt.scatter(merged['depth'], merged['median_loess_resid'], alpha=0.7)
    plt.axhline(0, color='grey', linestyle='--')
    plt.xlabel('Total Depth')
    plt.ylabel('Median LOESS Residual')
    plt.title(f'LOESS Residual Median vs. Depth (r = {r_loess:.2f})')
    plt.tight_layout()
    plt.savefig(
        "/Users/alexandra/PhD/PyCharmProjects/ALFAssay/DataSummary/plots/heatmap_bins/"+process+"_cov_ratio_loess.pdf")
    plt.clf()

def corr_depth_ratios(bins_df):
    depth_df = pd.read_csv("/Users/alexandra/PhD/healthy_sWGS/mosdepth/total_depth_healthy.txt", sep=" ", header=None)

    depth_df = depth_df.copy()
    depth_df.columns = ['raw_id', 'depth']
    depth_df['PatientId'] = depth_df['raw_id'].str.replace(
        r'\.mosdepth\.summary\.txt$', '', regex=True
    )
    depth_df = depth_df[['PatientId', 'depth']]

    # 2. compute per‐patient median ratio across the 204 bins
    #    (excluding the PatientId column)
    # corr_df = bins_df.merge(depth_df, on='PatientId', how='inner')
    #
    # bin_cols = [f'bin{i}' for i in range(1, 205)]
    #
    # # 4) Divide each bin by depth
    # corr_df[bin_cols] = corr_df[bin_cols].div(corr_df['depth'], axis=0)
    bins_df_log = bins_df.copy()
    for col in bins_df_log.columns.drop('PatientId'):
        bins_df_log[col] = np.log2(bins_df_log[col] + 1e-6)  # pseudocount

    median_ratio = bins_df_log.drop(columns='PatientId') \
        .median(axis=1)

    med_df = pd.DataFrame({
        'PatientId': bins_df_log['PatientId'],
        'median_ratio': median_ratio
    })

    # 3. merge with depth
    corr_df = med_df.merge(depth_df, on='PatientId', how='inner')

    # 4. compute Pearson r
    r_value = corr_df['median_ratio'].corr(corr_df['depth'])
    print(f'Pearson r = {r_value:.3f}')

    # optionally, get slope/intercept for plotting line
    slope, intercept, r, p, stderr = linregress(
        corr_df['depth'], corr_df['median_ratio']
    )

    # 5. plot
    plt.figure(figsize=(6,6))
    plt.scatter(corr_df['depth'], corr_df['median_ratio'])
    plt.plot(
        corr_df['depth'],
        intercept + slope * corr_df['depth'],
        linewidth=1
    )
    plt.xlabel('Total Depth')
    plt.ylabel('Median Ratio per Bin')
    plt.title(f'Median Ratio vs. Depth (r = {r_value:.2f})')
    plt.tight_layout()
    plt.savefig(
        "/Users/alexandra/PhD/PyCharmProjects/ALFAssay/DataSummary/plots/heatmap_bins/healthy_cov_ratio_log.pdf")
    plt.clf()


def fix_zeros(data_sorted, loess_frac=0.3):
    """
    Replace zeros in the DataFrame with a pseudocount and smooth profiles using LOESS.

    Parameters:
    - data_sorted: pd.DataFrame with PatientId column and numeric bin columns
    - loess_frac: float, the fraction of data used when estimating each y-value in LOESS

    Returns:
    - pd.DataFrame with PatientId and loess-smoothed values for each bin
    """
    # 1. Set PatientId as index
    df = data_sorted.set_index('PatientId')

    # 2. Replace zeros with NaN to identify entries to impute
    df = df.replace(0, np.nan)

    # 3. Compute pseudocount epsilon = 0.5 * smallest non-zero value
    min_nonzero = df[df > 0].min().min()
    epsilon = min_nonzero * 0.5

    # 4. Add pseudocount to all entries (removes zeros)
    df_pc = df + epsilon

    # 5. Mask original zeros for interpolation
    df_masked = df_pc.mask(df == 0, np.nan)

    # 6. Interpolate linearly across bins
    df_interp = df_masked.interpolate(
        axis=1,
        method='linear',
        limit_direction='both'
    )

    # 7. Forward/backward fill remaining NaNs at ends
    df_filled = df_interp.apply(
        lambda row: row.fillna(method='ffill').fillna(method='bfill'),
        axis=1
    )

    # 8. Smooth each profile using LOESS
    #    We use the column positions as the x-axis for smoothing
    bin_positions = np.arange(len(df_filled.columns))
    def _loess_smooth(row):
        smoothed = lowess(endog=row.values,
                          exog=bin_positions,
                          frac=loess_frac,
                          return_sorted=False)
        return smoothed

    loess_smoothed = df_filled.apply(_loess_smooth, axis=1)

    # 9. Reconstruct DataFrame with original columns and reset index
    df_smoothed = pd.DataFrame(
        loess_smoothed.tolist(),
        index=df_filled.index,
        columns=df_filled.columns
    ).reset_index()

    return df_smoothed


def frag_paterns_by_batch(file, task, study_analysed, features, metafile, bins, by_ichor):
    metadata = pd.read_csv(metafile, sep="\t")

    if study_analysed == "Synergy":
        metadata = metadata[~metadata["PatientId"].str.contains("Genome-IJB")]

    data = pd.read_csv(file, sep="\t")

    data = data.fillna(0)

    metadata = metadata[metadata[task].notna()]
    metadata["ctDNADetected"] = metadata["ctDNADetected"].apply(float)
    metadata["Label"] = metadata["Label"].apply(float)
    metadata["ctDNA"] = "ctDNA Positive"
    metadata["ctDNA"] = np.where(((metadata["ctDNADetected"] == 0) & (metadata["Label"] == 1)), "ctDNA Negative",
                                 metadata["ctDNA"])
    metadata["ctDNA"] = np.where(((metadata["ctDNADetected"] == 0) & (metadata["Label"] == 0)), "Healthy",
                                 metadata["ctDNA"])

    metadata["VAF"] = metadata["VAF"].apply(float)
    metadata = metadata[metadata["study"].isin(["Healthy", study_analysed])].sort_values(
        ["Label", "ichorTF", "batch", "PatientId"])

    if by_ichor:
        metadata["ctDNA"] = np.where(((metadata["ichorTF"] < 0.03) & (metadata["Label"] == 1)), "ctDNA Negative",
                                 metadata["ctDNA"])
        metadata["ctDNA"] = np.where(((metadata["ichorTF"] >= 0.03) & (metadata["Label"] == 1)), "ctDNA Positive",
                                     metadata["ctDNA"])

    missing = ["Genome-IJB-HP-40-xx_S28", "NIPT-PIJB-HP-73-Normal-15052023_S143"]
    meta_v2 = metadata[~metadata["PatientId"].isin(missing)]
    patients = list(meta_v2["PatientId"].values)
    patients_healthy = list(meta_v2[meta_v2["study"].isin(["Healthy"])]["PatientId"].values)
    patients_study = list(meta_v2[meta_v2["study"].isin([study_analysed])]["PatientId"].values)

    eps = 0.001

    data_feat_reads = data[data["Feature"] == "noReadsGCParagonFrag"][patients].apply(pd.to_numeric) + eps
    data_feat_short_reads = data[data["Feature"] == "shortFragGCPara"][patients].apply(pd.to_numeric)
    data_feat_ultra_long_reads = data[data["Feature"] == "ultraLongFragGCPara"][patients].apply(pd.to_numeric)

    ratio_long = data_feat_ultra_long_reads.values / data_feat_reads.values
    ratio_short = data_feat_short_reads.values / data_feat_reads.values

    ratio_healthy_long = data_feat_ultra_long_reads[patients_healthy].values / data_feat_reads[patients_healthy].values
    median_healthy_long = np.median(ratio_healthy_long, axis=1) + eps

    ratio_healthy_short = data_feat_short_reads[patients_healthy].values / data_feat_reads[patients_healthy].values
    median_healthy_short = np.median(ratio_healthy_short, axis=1) + eps

    ratio_study_long = data_feat_ultra_long_reads[patients_study].values / data_feat_reads[patients_study].values
    ratio_study_long_m = np.median(ratio_study_long, axis=1)

    ratio_study_short = data_feat_short_reads[patients_study].values / data_feat_reads[patients_study].values
    ratio_study_short_m = np.median(ratio_study_short, axis=1)

    cols_plot = ["bin" + str(bin) for bin in list(range(1, bins, 1))]

    for feature in features:

        # long
        t_test_study_healthy_long = kstest(median_healthy_long, ratio_study_long_m)
        # ratio_scaled_long = (ratio_long.T / median_healthy_long).T
        ratio_scaled_long = ratio_long
        data_sorted_long = pd.DataFrame(data=ratio_scaled_long.T, columns=cols_plot)
        data_sorted_long["PatientId"] = list(meta_v2["PatientId"].values)
        # data_sorted_long = fix_zeros(data_sorted_long)
        data_sorted_long[feature] = list(meta_v2[feature].values)
        data_sorted_long = data_sorted_long.sort_values([feature, "PatientId"])
        data_sorted_long = data_sorted_long.set_index("PatientId")

        labels = data_sorted_long.pop(feature)
        pallet = sns.cubehelix_palette(labels.unique().size,
                                       light=.9, dark=.1, reverse=True,
                                       start=1, rot=-2)
        lut = dict(zip(labels.unique(), pallet))
        row_colors = labels.map(lut)

        # labels_colors = pd.Series(labels, index=data_sorted_long.columns, name='labels').map(lut)

        # figsize = (18, 18)
        # fig, ax = plt.subplots(figsize=figsize) Genome-IJB-HP-70_S175

        # sns.heatmap(ax=ax, data=data_sorted_long.T)
        plot_long = sns.clustermap(data_sorted_long[cols_plot], row_colors=row_colors, row_cluster=False,
                                   col_cluster=False)
        # plt.show()
        # plt.title("Long fragments on bins for healthy and " + study_analysed + ", t-test=" + str(round(t_test_study_healthy_long[1],2)))
        title = "Long fragments on bins for healthy and " + study_analysed + ", t-test=" + str(
            round(t_test_study_healthy_long[1], 2))
        u.fix_clustermap(plot_long, title)
        for label in labels.unique():
            plot_long.ax_col_dendrogram.bar(0, 0, color=lut[label], label=label, linewidth=0)
        l1 = plot_long.ax_col_dendrogram.legend(title='Status', loc="center", ncol=5, bbox_to_anchor=(0.35, 0.89),
                                                bbox_transform=gcf().transFigure)

        # plt.show()
        plt.savefig(
            "/Users/alexandra/PhD/PyCharmProjects/ALFAssay/DataSummary/plots/heatmap_bins/long_healthy_" + study_analysed + "_on_" + feature + ".pdf")
        plt.clf()

        ##short
        t_test_study_healthy_short = kstest(median_healthy_short, ratio_study_short_m)
        # ratio_scaled_short = (ratio_short.T / median_healthy_short).T
        ratio_scaled_short = ratio_short
        data_sorted_short = pd.DataFrame(data=ratio_scaled_short.T, columns=cols_plot)
        data_sorted_short["PatientId"] = list(meta_v2["PatientId"].values)
        # data_sorted_short = fix_zeros(data_sorted_short)
        data_sorted_short[feature] = list(meta_v2[feature].values)
        data_sorted_short = data_sorted_short.sort_values([feature, "PatientId"])
        data_sorted_short = data_sorted_short.set_index("PatientId")
        # figsize = (18, 18)
        # fig, ax = plt.subplots(figsize=figsize)

        # sns.heatmap(ax=ax, data=data_sorted_short.T)

        plot_short = sns.clustermap(data_sorted_short[cols_plot], row_colors=row_colors, row_cluster=False,
                                    col_cluster=False)
        title = "Short fragments on bins for healthy and " + study_analysed + ", t-test=" + str(
            round(t_test_study_healthy_short[1], 2))
        u.fix_clustermap(plot_short, title)
        for label in labels.unique():
            plot_short.ax_col_dendrogram.bar(0, 0, color=lut[label], label=label, linewidth=0)
        l1 = plot_short.ax_col_dendrogram.legend(title='Status', loc="center", ncol=3, bbox_to_anchor=(0.35, 0.89),
                                                 bbox_transform=gcf().transFigure)

        plt.savefig(
            "/Users/alexandra/PhD/PyCharmProjects/ALFAssay/DataSummary/plots/heatmap_bins/short_healthy_" + study_analysed + "_on_" + feature + ".pdf")
        plt.clf()


# studies = ["Pearl", "Synergy", "Neorhea"]
studies = ["Synergy"]
file = "/Users/alexandra/PhD/FragmentationPatterns/Data/NN_5MB/NN_5MB_pivot.csv"
task = "ctDNADetected"
# meta_file = "/Users/alexandra/PhD/PyCharmProjects/ALFAssayNN/results/results_validation_PRL.csv"
meta_file = "/Users/alexandra/PhD/FragmentationPatterns/Data/MetaData/AllMetaData.csv"
# features = ["Label", "batch", "ichorTF", "study", "ctDNADetected"]
# features = ["Predictedlabel"]
features = ["ctDNA"]
bins = 205
by_ichor = True
# for study in studies:
#     frag_paterns_by_batch(file, task, study, features, meta_file, bins, by_ichor)



def frag_paterns_bins(file, task, features, metafile, bins, by_ichor,process):
    metadata_ = pd.read_csv(metafile, sep="\t")

    data = pd.read_csv(file, sep="\t")

    data = data.fillna(0)

    metadata = metadata_[metadata_[task].notna()]
    # if process=="Training":
    #     lowcov = ["NIPT-P0104-ArmB-SNRFR0085-W1-Screen_S37", "NIPT-P0154-ArmA-SNRFR0005-W3_S44",
    #               "NIPT-P0156-ArmA-SNRFR0046-W1-Screen_S46", "NIPT-P0156-ArmA-SNRFR0046-W13-EOC_S48",
    #               "NIPT-P0156-ArmA-SNRFR0046-W3_S47"]  #
    #
    #     # metadata = metadata[(metadata["study"]!="Neorhea") ]
    #     metadata = metadata[~metadata["PatientId"].isin(lowcov)]
    metadata = metadata[metadata["Process"]==process]
    metadata["ctDNADetected"] = metadata["ctDNADetected"].apply(float)
    metadata["Label"] = metadata["Label"].apply(float)
    metadata["ctDNA"] = "ctDNA Positive"
    metadata["ctDNA"] = np.where(((metadata["ctDNADetected"] == 0) & (metadata["Label"] == 1)), "ctDNA Negative",
                                 metadata["ctDNA"])
    metadata["ctDNA"] = np.where(((metadata["ctDNADetected"] == 0) & (metadata["Label"] == 0)), "Healthy",
                                 metadata["ctDNA"])

    metadata["VAF"] = metadata["VAF"].apply(float)
    metadata = metadata.sort_values( ["Label", "ichorTF", "batch", "PatientId"])

    if by_ichor:
        metadata["ctDNA"] = np.where(((metadata["ichorTF"] < 0.03) & (metadata["Label"] == 1)), "ctDNA Negative",
                                 metadata["ctDNA"])
        metadata["ctDNA"] = np.where(((metadata["ichorTF"] >= 0.03) & (metadata["Label"] == 1)), "ctDNA Positive",
                                     metadata["ctDNA"])

    # missing = ["Genome-IJB-HP-40-xx_S28", "NIPT-PIJB-HP-73-Normal-15052023_S143"]
    # meta_v2 = metadata[~metadata["PatientId"].isin(missing)]
    patients = list(metadata[metadata["Process"]==process]["PatientId"].values)
    # patients_healthy = list(meta_v2[meta_v2["study"].isin(["Healthy"])]["PatientId"].values)
    # patients_study = list(meta_v2[meta_v2["study"].isin([study_analysed])]["PatientId"].values)


    ratio_long = data[data["Feature"] == "ratioLong"][patients].apply(pd.to_numeric)
    ratio_short = data[data["Feature"] == "ratioShort"][patients].apply(pd.to_numeric)


    cols_plot = ["bin" + str(bin) for bin in list(range(1, bins, 1))]
    cols_short = ["bin_short_" + str(bin) for bin in list(range(1, bins, 1))]
    cols_long = ["bin_long_" + str(bin) for bin in list(range(1, bins, 1))]


    for feature in features:

        # long
        # t_test_study_healthy_long = kstest(median_healthy_long, ratio_study_long_m)
        # ratio_scaled_long = (ratio_long.T / median_healthy_long).T
        # ratio_scaled_long = ratio_long
        # data_sorted_long = pd.DataFrame(data=ratio_scaled_long.T, columns=cols_plot)
        # data_sorted_long["PatientId"] = list(metadata["PatientId"].values)
        # # data_sorted_long = fix_zeros(data_sorted_long)
        # data_sorted_long[feature] = list(metadata[feature].values)

        df_long = ratio_long.T.reset_index()

        df_long.columns = ["PatientId"] + cols_plot
        # df_long = fix_zeros(df_long)
        df_long = df_long.merge(metadata[["PatientId"]+[feature]], on="PatientId", how="inner")


        data_sorted_long = df_long.sort_values([feature, "PatientId"])
        data_sorted_long = data_sorted_long.set_index("PatientId")

        # df = ratio_long.T.reset_index()
        # df.columns = ["PatientId"] + cols_plot
        # df_long = df.merge(metadata_, on="PatientId", how="left")

        labels = data_sorted_long.pop(feature)
        pallet = sns.cubehelix_palette(labels.unique().size,
                                       light=.9, dark=.1, reverse=True,
                                       start=1, rot=-2)
        custom_colors = ["#56B4E9", "#009E73", "#E69F00"]
        lut = dict(zip(labels.unique(), custom_colors))
        row_colors = labels.map(lut)
        cmap = sns.diverging_palette(250, 10, as_cmap=True)
        # labels_colors = pd.Series(labels, index=data_sorted_long.columns, name='labels').map(lut)

        # figsize = (18, 18)
        # fig, ax = plt.subplots(figsize=figsize) Genome-IJB-HP-70_S175

        # sns.heatmap(ax=ax, data=data_sorted_long.T)
        plot_long = sns.clustermap(data_sorted_long[cols_plot], row_colors=row_colors, row_cluster=False,
                                   col_cluster=False, cmap=cmap)
        # plt.show()
        # plt.title("Long fragments on bins for healthy and " + study_analysed + ", t-test=" + str(round(t_test_study_healthy_long[1],2)))
        title = "Long fragments on bins for  " + process + " dataset"
        u.fix_clustermap(plot_long, title)
        for label in labels.unique():
            plot_long.ax_col_dendrogram.bar(0, 0, color=lut[label], label=label, linewidth=0)
        l1 = plot_long.ax_col_dendrogram.legend(title='Status', loc="center", ncol=5, bbox_to_anchor=(0.35, 0.89),
                                                bbox_transform=gcf().transFigure)

        # plt.show()
        plt.savefig(
            "/Users/alexandra/PhD/PyCharmProjects/ALFAssay/DataSummary/plots/heatmap_bins/long_ratio_on_" + feature + "_" + process + "_loes.pdf")
        plt.clf()

        ##short
        # t_test_study_healthy_short = kstest(median_healthy_short, ratio_study_short_m)
        # ratio_scaled_short = (ratio_short.T / median_healthy_short).T
        # ratio_scaled_short = ratio_short
        # data_sorted_short = pd.DataFrame(data=ratio_scaled_short.T, columns=cols_plot)
        # data_sorted_short["PatientId"] = list(metadata["PatientId"].values)
        # # data_sorted_short.to_csv("data_sorted_short_training.csv", index=False, sep="\t")
        # # data_sorted_short = fix_zeros(data_sorted_short)
        # data_sorted_short[feature] = list(metadata[feature].values)

        df_short = ratio_short.T.reset_index()

        df_short.columns = ["PatientId"] + cols_plot
        # if process=="Training":
        corr_depth_loess(df_short, metadata, process)
        # df_short = fix_zeros(df_short)
        df_short = df_short.merge(metadata[["PatientId","ichorTF"] + [feature]], on="PatientId", how="inner")


        data_sorted_short = df_short.sort_values([feature, "PatientId"])
        data_sorted_short = data_sorted_short.set_index("PatientId")

        # df = data[data["Feature"] == "ratioShort"].T.reset_index()
        # df.columns = ["PatientId"] + cols_plot
        # df_short = df.merge(metadata_, on ="PatientId", how="left")
        #
        df_short.to_csv("data_sorted_short_"+process+".csv", index=False, sep="\t")

        # figsize = (18, 18)
        # fig, ax = plt.subplots(figsize=figsize)

        # sns.heatmap(ax=ax, data=data_sorted_short.T)

        plot_short = sns.clustermap(data_sorted_short[cols_plot], row_colors=row_colors, row_cluster=False,
                                    col_cluster=False,cmap=cmap)
        title = "Short fragments on bins for healthy and " + str(process) + " dataset"
        u.fix_clustermap(plot_short, title)
        for label in labels.unique():
            plot_short.ax_col_dendrogram.bar(0, 0, color=lut[label], label=label, linewidth=0)
        l1 = plot_short.ax_col_dendrogram.legend(title='Status', loc="center", ncol=3, bbox_to_anchor=(0.35, 0.89),
                                                 bbox_transform=gcf().transFigure)

        plt.savefig(
            "/Users/alexandra/PhD/PyCharmProjects/ALFAssay/DataSummary/plots/heatmap_bins/short_healthy_on_" + feature + "_"+ process + "_loes.pdf")
        plt.clf()


# studies = ["Pearl", "Synergy", "Neorhea"]
processes = ["Training", "Validation"]
file = "/Users/alexandra/PhD/FragmentationPatterns/Data/NN_5MB/NN_5MB_pivot_with_ratios_.csv"
task = "ctDNADetected"
# meta_file = "/Users/alexandra/PhD/PyCharmProjects/ALFAssayNN/results/results_validation_PRL.csv"
meta_file = "/Users/alexandra/PhD/FragmentationPatterns/Data/MetaData/AllMetaData_.csv"
# features = ["Label", "batch", "ichorTF", "study", "ctDNADetected"]
# features = ["Predictedlabel"]
features = ["ctDNA"]
bins = 205
by_ichor = True
for process in processes:
    frag_paterns_bins(file, task, features, meta_file, bins, by_ichor, process)
