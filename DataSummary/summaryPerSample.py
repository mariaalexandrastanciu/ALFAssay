# Created by alexandra at 02/08/2023
import pandas as pd
from FeatureEngineering import utils as u
import glob
import os
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats as stats
from scipy.stats import ks_2samp

def calc_summaries_per_sample(file, meta_data_file, name_pattern):
    patient_id = u.get_patient_from_file_name(file, name_pattern)
    df = pd.read_csv(file, sep="\t")
    meta_data = pd.read_csv(meta_data_file, sep="\t")
    meta_data_sample = meta_data[meta_data["PatientId"] == patient_id].values.tolist()[0]

    median_frag_analysed = int(df["fragment_size"].median())
    max_count_frag_analysed = int(df[df["fragment_counts"] == df["fragment_counts"].max()]["fragment_size"].values[0])
    mean_frag_analysed = int(df["fragment_size"].mean())
    q25_frag_analysed = int(df["fragment_size"].quantile(.25))
    q75_frag_analysed = int(df["fragment_size"].quantile(.75))
    total_frag_analysed = int(df["fragment_counts"].sum())

    summary = ["Median", "MaxCountSize", "Mean", "q25", "q75", "TotalNo"]
    cols = summary + meta_data.columns.values.tolist()
    summary_data = [median_frag_analysed,
                    max_count_frag_analysed,
                    mean_frag_analysed,
                    q25_frag_analysed,
                    q75_frag_analysed,
                    total_frag_analysed] + meta_data_sample
    summary_df = pd.DataFrame(data=[summary_data],
                              columns=cols)

    ## delfi comparison?
    delfi_frags = df[(df["fragment_size"] >= 100) & (df["fragment_size"] <= 220)]
    median_delfi_frag = int(delfi_frags["fragment_counts"].median())
    max_count_frag_delfi_frag = int(delfi_frags[delfi_frags["fragment_counts"] ==
                                                delfi_frags["fragment_counts"].max()]["fragment_size"].values[0])
    mean_delfi_frag = int(delfi_frags["fragment_size"].mean())
    q25_delfi_frag = int(delfi_frags["fragment_size"].quantile(.25))
    q75_delfi_frag = int(delfi_frags["fragment_size"].quantile(.75))
    total_delfi_frags = int(delfi_frags["fragment_counts"].sum())

    summary_delfi_data = [median_delfi_frag,
                          max_count_frag_delfi_frag,
                          mean_delfi_frag,
                          q25_delfi_frag,
                          q75_delfi_frag,
                          total_delfi_frags] + meta_data_sample
    summary_delfi_df = pd.DataFrame(data=[summary_delfi_data],
                                    columns=cols)

    return summary_df, summary_delfi_df


def calc_summary_by_study(path_to_files, meta_data_file, output_file, name_pattern):
    frag_size_files = glob.glob(path_to_files + "/*.csv")
    summary_df = pd.DataFrame()
    summary_delfi_df = pd.DataFrame()
    for frag_size_file in frag_size_files:
        sample_summary, sample_summary_delfi = calc_summaries_per_sample(frag_size_file, meta_data_file, name_pattern)
        summary_df = pd.concat([summary_df, sample_summary], ignore_index=True)
        summary_delfi_df = pd.concat([summary_delfi_df, sample_summary_delfi], ignore_index=True)
    summary_df.to_csv(output_file, sep="\t")
    output_file_delfi = output_file.replace(".csv", "") + "_delfi.csv"
    summary_delfi_df.to_csv(output_file_delfi, sep="\t")


def plot_density(file, plot_output, start_frag_size, end_frag_size, title):
    cols_with_label = [str(i) for i in range(start_frag_size, end_frag_size)] + [
        "Label"]  # list(range(0, nr_of_frag_sizes))
    cols = [str(i) for i in range(start_frag_size, end_frag_size)]
    sns.set_style('whitegrid')
    dataset = pd.read_csv(file, sep="\t")
    if "Synergy" in title:
        dataset = dataset[~dataset["PatientId"].str.contains("Genome-IJB")]
    data = dataset[cols_with_label].groupby(["Label"])[cols].median().T
    data_scaled = data.div(data.sum(axis=0), axis=1)  # minmax_scale(data, feature_range=(0, 1), axis=0, copy=True)
    data_scaled = pd.DataFrame(data_scaled, columns=['Cancer', "Healthy"])
    frag_sizes = list(range(start_frag_size, end_frag_size))
    data_scaled['FragmentSizes'] = frag_sizes
    data_melted = pd.melt(data_scaled, id_vars=['FragmentSizes'], value_vars=['Cancer', "Healthy"])
    data_melted.rename(columns={"value": "FragmentSizeCounts", "variable": "Type"}, inplace=True)

    sns.lineplot(x="FragmentSizes", y="FragmentSizeCounts", hue="Type", data=data_melted)
    plt.xlabel("Fragment Sizes")
    plt.title(title, fontsize=8)
    plt.savefig(plot_output)
    plt.clf()


def calc_summary(file_raw, start_frag_size, end_frag_size, output_file):
    dataset_raw = pd.read_csv(file_raw, sep="\t")
    cols_with_label = [str(i) for i in range(start_frag_size, end_frag_size - start_frag_size)] + ["Label"]
    cols_delfi = [str(i) for i in range(100, 220)] + ["Label"]
    data_raw = dataset_raw[cols_with_label]
    delfi_data_raw = data_raw[cols_delfi]
    data_raw['TotalNrFragment'] = data_raw.sum(axis=1, numeric_only=True)
    delfi_data_raw['TotalNrFragmentDELFI'] = delfi_data_raw.sum(axis=1, numeric_only=True)

    summary_cols = ["PatientId", "Label", "ichorTF", "VAF", "ctDNADetected", "survivalStatus", "timepoint",
                    "TotalNrFragmentRaw", "TotalNrFragmentDELFIRaw"]
    summary_df = pd.DataFrame([dataset_raw["PatientId"].values, dataset_raw["Label"].values,
                               dataset_raw["ichorTF"].values, dataset_raw["VAF"].values,
                               dataset_raw["ctDNADetected"].values, dataset_raw["survivalStatus"].values,
                               dataset_raw["timepoint"].values, data_raw['TotalNrFragment'].values,
                               delfi_data_raw['TotalNrFragmentDELFI'].values
                               ]).T
    summary_df.columns = summary_cols

    ## columns with frag sizes ranges ##
    cols_very_short_frag = [str(i) for i in range(start_frag_size, 100)]  ## pos 30 - 100
    cols_short_frag = [str(i) for i in range(90 , 150  )]  ## pos 100 - 150
    cols_long_frag = [str(i) for i in range(151  , 220  )]  ## pos 151 - 220
    cols_short_dinu_frag = [str(i) for i in range(221  , 300 )]  ## 221 - 300
    cols_long_dinu_frag = [str(i) for i in range(301  , 360 )]  ## 301 - 360
    cols_short_trinu_frag = [str(i) for i in range(361 , 500  )]  ## 361 -500
    cols_long_trinu_frag = [str(i) for i in range(501 , 560  )]  ## 501 - 561
    cols_very_long_frag = [str(i) for i in range(561  , 700  )]  ## 561 - 700
    cols_density_plot = [str(i) for i in range(300, 400)]  ## 561 - 700
    cols_long = [str(i) for i in range(300, 700)]

    cols_frag_sizes = [str(i) for i in range(start_frag_size, end_frag_size )]
    ## raw data ##
    dataset_raw_scaled = dataset_raw[cols_frag_sizes].div(dataset_raw[cols_frag_sizes].sum(axis=1), axis=0)
    summary_df["very_short_30_100_Raw"] = dataset_raw_scaled[cols_very_short_frag].sum(axis=1)
    summary_df["short_n_90_150_Raw"] = dataset_raw_scaled[cols_short_frag].sum(axis=1)
    summary_df["long_n_151_220_Raw"] = dataset_raw_scaled[cols_long_frag].sum(axis=1)
    summary_df["short_din_221_300_Raw"] = dataset_raw_scaled[cols_short_dinu_frag].sum(axis=1)
    summary_df["long_din_301_360_Raw"] = dataset_raw_scaled[cols_long_dinu_frag].sum(axis=1)
    summary_df["short_trin_361_500_Raw"] = dataset_raw_scaled[cols_short_trinu_frag].sum(axis=1)
    summary_df["long_trin_501_561_Raw"] = dataset_raw_scaled[cols_long_trinu_frag].sum(axis=1)
    summary_df["very_long_561_700_Raw"] = dataset_raw_scaled[cols_very_long_frag].sum(axis=1)
    summary_df["long_300_400_Raw"] = dataset_raw_scaled[cols_density_plot].sum(axis=1)
    summary_df["long_300_700"] = dataset_raw_scaled[cols_long].sum(axis=1)

    summary_df.to_csv(output_file, sep="\t")


def plot_boxplot(summary_file, output_file, title):
    dataset = pd.read_csv(summary_file, sep="\t")
    sns.set(style="darkgrid", color_codes=True)
    dataset_boxplot = dataset[["PatientId", "TotalNrFragmentRaw", "TotalNrFragmentGCCorrected", "Label"]]
    data_melted = pd.melt(dataset_boxplot, id_vars=['PatientId', "Label"], value_vars=["TotalNrFragmentRaw",
                                                                                       "TotalNrFragmentGCCorrected"])
    data_melted = data_melted.rename(columns={"value": "No of fragments"})
    sns.boxplot(data=data_melted, x="Label", y="No of fragments", hue="variable")
    plt.title(title, fontsize=6)
    plt.savefig(output_file)
    plt.clf()


# def plot_boxplot_between_studies(summary_file_1, summary_file_2, output_file, title, studies):
#     dataset_1 = pd.read_csv(summary_file_1, sep="\t")
#     dataset_2 = pd.read_csv(summary_file_2, sep="\t")
#     sns.set(style="darkgrid", color_codes=True)
#     dataset_1 =
#     dataset = pd.DataFrame([studies + ])
#     sns.boxplot(data=data_melted, x="Label", y="value", hue="variable")


def plot_response_corr_with_frag_size(file, output_file, cols, title):
    dataset = pd.read_csv(file, sep="\t")
    dataset_cancer = dataset[dataset["Label"] == "Cancer"]
    for col in cols:
        fig, axes = plt.subplots(2, 3, figsize=(18, 10))
        print("Create plot for frag sizes:" + col)
        fig.suptitle("Correlation for fragments sizes: " + col + ";" + title)
        ## Label ##
        label_test = stats.mannwhitneyu(dataset[dataset["Label"] == "Cancer"][col],
                                     dataset[dataset["Label"] == "Healthy"][col])
        sns.boxplot(ax=axes[0, 0], data=dataset, x=col, y='Label')
        axes[0, 0].title.set_text("Mann Whitney U Test p=" + str(np.round(label_test.pvalue, 4)))

        ## ctDNA detection ##
        ## some bad thing: TODO: rethink this ##
        if "Pearl" in file:
            dataset.loc[dataset["VAF"].isna(), "ctDNADetected"] = False
            dataset.loc[dataset["VAF"] > 0, "ctDNADetected"] = True
        ###
        dataset_ctDNADetection = dataset[["Label", "ctDNADetected", col]].dropna()
        dataset_ctDNADetection.loc[dataset_ctDNADetection["ctDNADetected"] == True, "ctDNADetected"] = "ctDNA Detected"
        dataset_ctDNADetection.loc[dataset_ctDNADetection["ctDNADetected"] == False, "ctDNADetected"] = \
            "ctDNA not Detected"
        dataset_ctDNADetection_cancer = dataset_ctDNADetection[dataset_ctDNADetection["Label"]=="Cancer"]
        if len(dataset_ctDNADetection_cancer[
                   dataset_ctDNADetection_cancer["ctDNADetected"] == "ctDNA Detected"])>0\
                and len(dataset_ctDNADetection_cancer[
                            dataset_ctDNADetection_cancer["ctDNADetected"] == "ctDNA not Detected"])>0:

            ctDNADetection_test = stats.mannwhitneyu(
                dataset_ctDNADetection_cancer[dataset_ctDNADetection_cancer["ctDNADetected"] == "ctDNA Detected"][col],
                dataset_ctDNADetection_cancer[dataset_ctDNADetection_cancer["ctDNADetected"] == "ctDNA not Detected"][col])
            ctDNADetection_test = ctDNADetection_test.pvalue
        else:
            ctDNADetection_test = 1

        sns.boxplot(ax=axes[0, 1], data=dataset_ctDNADetection, x=col, y='ctDNADetected', hue="Label")
        axes[0, 1].title.set_text("Mann Whitney U Test  p=" + str(np.round(ctDNADetection_test, 4)))

        ## survival test ###
        dataset_survivalStatus = dataset[["Label", "survivalStatus", col]].dropna()
        dataset_survivalStatus.loc[dataset_survivalStatus["survivalStatus"] == 1, "survivalStatus"] = "Relapse"
        dataset_survivalStatus.loc[dataset_survivalStatus["survivalStatus"] == 0, "survivalStatus"] = "No relapse"
        dataset_survivalStatus_cancer = dataset_survivalStatus[dataset_survivalStatus["Label"] == "Cancer"]

        if len(dataset_survivalStatus_cancer[
                   dataset_survivalStatus_cancer["survivalStatus"] == "Relapse"]) > 0:
            survival_test = stats.mannwhitneyu(
                dataset_survivalStatus_cancer[dataset_survivalStatus_cancer["survivalStatus"] == "Relapse"][col],
                dataset_survivalStatus_cancer[dataset_survivalStatus_cancer["survivalStatus"] == "No relapse"][col])
            survival_test = survival_test.pvalue
        else:
            survival_test = 1
        sns.boxplot(ax=axes[0, 2], data=dataset_survivalStatus, x=col, y='survivalStatus', hue="Label")
        axes[0, 2].title.set_text("Mann Whitney U Test  p=" + str(np.round(survival_test, 4)))
        #### ichorTF ##
        ichortf_corr = stats.pearsonr(dataset_cancer[col], dataset_cancer["ichorTF"])
        sns.scatterplot(ax=axes[1, 0], data=dataset, x=col, y='ichorTF', hue="Label")
        axes[1, 0].title.set_text("Correlation p=" + str(np.round(ichortf_corr[1], 4)))
        ### VAF ##

        dataset_VAF = dataset[["Label", "VAF", col]].dropna()
        dataset_VAF_cancer = dataset_cancer[["Label", "VAF", col]].dropna()
        if dataset_VAF_cancer.empty:
            ichortf_corr = stats.pearsonr(dataset_cancer[col], dataset_cancer["ichorTF"])
            sns.scatterplot(ax=axes[1, 0], data=dataset, x=col, y='ichorTF', hue="Label")
            axes[1, 1].title.set_text("Correlation p=" + str(np.round(ichortf_corr[1], 4)))
        else:
            vaf_corr = stats.pearsonr(dataset_VAF_cancer[col], dataset_VAF_cancer["VAF"])
            sns.scatterplot(ax=axes[1, 1], data=dataset_VAF, x=col, y='VAF', hue="Label")
            axes[1, 1].title.set_text("Correlation p=" + str(np.round(vaf_corr[1], 4)))

        ## timepoint ###
        timepoint_test = stats.kruskal(dataset_cancer[dataset_cancer["timepoint"] == "BL"][col],
                                       dataset_cancer[dataset_cancer["timepoint"] == "C1"][col],
                                       dataset_cancer[dataset_cancer["timepoint"] == "Surgery"][col])
        sns.boxplot(ax=axes[1, 2], data=dataset_cancer, x=col, y='timepoint')
        axes[1, 2].title.set_text("kruskal wallis test p=" + str(np.round(timepoint_test.pvalue, 4)))
        # for ax in axes:
        #     for ax_line in ax:
        #         # ax_line.text(0.85, 0.85, ctDNADetection_test[1], fontsize=9)

        plt.subplots_adjust(bottom=0.2, right=0.8, top=0.9, wspace=0.8)
        save_file = output_file + "_" + col + ".pdf"
        plt.savefig(save_file)
        plt.clf()

def mutations_summary(path_to_mutations_plasma, path_to_mutation_tumor, output_path):
    files = glob.glob(path_to_mutations_plasma + "/*.csv")
    cols = ["SampleId", "MutationNrPlasma", "MutationNrTumor"]
    mutation_summary = []
    for file in files:
        file_name = os.path.basename(file).replace(".csv", "")
        data = pd.read_csv(file, sep="\t")
        data_tumor = pd.read_csv(path_to_mutation_tumor + "/" + file_name + ".csv", sep="\t")
        mutation_summary.append([file_name, len(data.index), len(data_tumor.index)])
    mutation_summary_df = pd.DataFrame(mutation_summary, columns=cols)
    save_file = output_path + "/summary_nr_mutations.csv"
    mutation_summary_df.to_csv(save_file, index=False, sep="\t")

    sns.histplot(data=mutation_summary_df, x="MutationNrPlasma")
    save_plot = output_path + "/Plots/nr_mutation_plasma_histogram.pdf"
    plt.savefig(save_plot)
    plt.clf()

    sns.histplot(data=mutation_summary_df, x="MutationNrTumor")
    save_plot = output_path + "/Plots/nr_mutation_tumor_histogram.pdf"
    plt.savefig(save_plot)
    plt.clf()

    mutation_summary_melt = pd.melt(mutation_summary_df, id_vars=['SampleId'], value_vars=['MutationNrPlasma',
                                                                                           "MutationNrTumor"])
    mutation_summary_melt.rename(columns={"value": "No mutations", "variable": "Mutation Source"}, inplace=True)

    ax = sns.boxplot(x='Mutation Source', y='No mutations', data=mutation_summary_melt)

    patients_with_many_mutations = mutation_summary_melt[(mutation_summary_melt["Mutation Source"] == "MutationNrPlasma")
                                                         & (mutation_summary_melt["No mutations"] > 50)]["SampleId"]

    mutation_summary_melt_filt = mutation_summary_melt[mutation_summary_melt["SampleId"].
                                                       isin(patients_with_many_mutations)]

    ax = sns.swarmplot(x='Mutation Source', y='No mutations', data=mutation_summary_melt_filt)
    sns.lineplot(data=mutation_summary_melt_filt, x="Mutation Source", y="No mutations", units="SampleId",
                 estimator=None, color=".5")
    label_test = stats.mannwhitneyu(mutation_summary_melt[mutation_summary_melt["Mutation Source"] == "MutationNrPlasma"]
                                    ["No mutations"],
                                    mutation_summary_melt[mutation_summary_melt["Mutation Source"] == "MutationNrTumor"]
                                    ["No mutations"])
    ax.title.set_text("Mann Whitney U Test p=" + str(np.round(label_test.pvalue, 4)))

    plt.yscale('log')
    save_plot = output_path + "/Plots/mutation_boxplot.pdf"
    plt.savefig(save_plot)
    plt.clf()


# create summary #########
def create_summary(input_files, output_summary_dir):
    start_frag_size_v = 30
    end_frag_size_v = 700
    raw_files = glob.glob(input_files + "/*")
    for file in raw_files:
        raw_file_name = os.path.basename(file).replace(".csv", "")
        print("Run analysis for: " + raw_file_name)
        summary_output_file = output_summary_dir + raw_file_name + ".csv"
        calc_summary(file, start_frag_size_v, end_frag_size_v, summary_output_file)


# plot fragment density ###
def plot_fragment_density(density_plot_output_dir, input_files):
    start_frag_size = 30
    end_frag_size = 700
    files = glob.glob(input_files)

    for file in files:
        file_name = os.path.basename(file).replace(".csv", "")
        print("Run analysis for: " + file_name)
        density_plot_output = density_plot_output_dir + file_name + ".pdf"
        title = "Fragment sizes distribution: " + u.get_title_from_filename(file_name)
        plot_density(file, density_plot_output,  start_frag_size, end_frag_size,  title)


# plot total number of fragments ###
def plot_total_nr_of_frag(output_dir_plot, input_files):
    files = glob.glob(input_files)
    for file in files:
        file_name = os.path.basename(file).replace(".csv", "")
        print("Run analysis for: " + file_name)
        corr_plot_output = output_dir_plot + file_name + ".pdf"
        title = "Absolute value of no of fragments for " + u.get_title_from_filename(file_name) + " compared with healthy"
        plot_boxplot(file, corr_plot_output, title)


# plot correlation ##
def plot_correlations(root):
    output_plot = root + "/Plots/FragSizesCountCorrelations/"
    files = glob.glob(root + "/*_fragCounts.csv")

    for file in files:
        file_name = os.path.basename(file).replace(".csv", "")
        if os.path.exists(root + "/Plots/FragSizesCountCorrelations/" + file_name):
            corr_plot_output = output_plot + file_name
        else:
            os.makedirs(root + "/Plots/FragSizesCountCorrelations/" + file_name)
            corr_plot_output = output_plot + file_name
        print("Run analysis for: " + file_name)
        # corr_plot_output = output_plot + file_name
        title = "File:" + u.get_title_from_filename(file_name)
        cols_raw = ["very_short_30_100_Raw", "short_n_90_150_Raw", "long_n_151_220_Raw", "short_din_221_300_Raw",
                "long_din_301_360_Raw", "short_trin_361_500_Raw", "long_trin_501_561_Raw", "very_long_561_700_Raw",
                "long_300_400_Raw", "long_300_700"]

        corr_plot_output_raw = corr_plot_output + file_name
        plot_response_corr_with_frag_size_2(file, corr_plot_output_raw, cols_raw, title)


def mutations_summary_test(path_to_mutations_plasma, output_path, title):
    files = glob.glob(path_to_mutations_plasma + "/*.csv")
    cols = ["SampleId", "MutationNrPlasma"]
    mutation_summary = []
    for file in files:
        file_name = os.path.basename(file).replace(".csv", "")
        data = pd.read_csv(file, sep="\t")
        mutation_summary.append([file_name, len(data.index)])
    mutation_summary_df = pd.DataFrame(mutation_summary, columns=cols)
    save_file = output_path + "/" + title + "summary_nr_mutations.csv"
    mutation_summary_df.to_csv(save_file, index=False, sep="\t")

    sns.histplot(data=mutation_summary_df, x="MutationNrPlasma")
    plt.title(title)
    save_plot = output_path + "/Plots/" + title + "nr_mutation_plasma_histogram.pdf"
    plt.savefig(save_plot)
    plt.clf()

def test_mutations_in_healthy(file1,file2):
    x = pd.read_csv(file1, sep="\t")["MutationNrPlasma"].values
    y = pd.read_csv(file2, sep="\t")["MutationNrPlasma"].values
    test = ks_2samp(x, y)
    print(test)
# mutations
def summary_histo_mutations():
    outpath_mutation_no = "/Users/alexandra/PhD/FragmentationPatterns/Summary/Mutations"
    path_to_mutations_per_sample = "/Users/alexandra/PhD/NeoRheaStudy/FragmentationPatterns/MutScanParssed"
    path_to_mutations_per_tumor = "/Users/alexandra/PhD/NeoRheaStudy/FragmentationPatterns/MutationsPerSample"
    mutations_summary(path_to_mutations_per_sample, path_to_mutations_per_tumor, outpath_mutation_no)


def plot_Neorh_prop_with_mutations():
    neorhea_data = pd.read_csv("/Users/alexandra/PhD/FragmentationPatterns/Summary/NeoRheaStudy_fragCounts.csv",
                               sep="\t")
    mutations_data = pd.read_csv("/Users/alexandra/PhD/FragmentationPatterns/Summary/Mutations/summary_nr_mutations.csv"
                                 , sep="\t")
    neorhea_data.rename(columns={"PatientId": "SampleId"}, inplace=True)
    all_data = pd.merge(neorhea_data, mutations_data, how="inner", on=["SampleId"])
    summary_cols_scatterplot_labels = ["ichorTF", "VAF", "TotalNrFragmentRaw", "TotalNrFragmentDELFIRaw",
                                       "TotalNrFragmentGCCorrected", "TotalNrFragmentDELFIGCCorrected"]
    summary_cols_scatterplot_frag_size = ["very_short_30_100_GC", "short_n_100_150_GC",
                                "long_n_151_220_GC", "short_din_221_300_GC", "long_din_301_360_GC",
                                "short_trin_361_500_GC", "long_trin_501_561_GC", "very_long_561_700_GC"]
    summary_cols_boxplot = ["ctDNADetected", "survivalStatus", "timepoint"]
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    i = 0
    j = 0
    print("Create plots for no of mutations in plasma and continuous labels")
    for col in summary_cols_scatterplot_labels:
        fig.suptitle(print("Create plots for no of mutations in plasma and continuous labels"))
        dataset_VAF = all_data[["MutationNrPlasma", col]].dropna()
        vaf_corr = stats.pearsonr(dataset_VAF[col], dataset_VAF["MutationNrPlasma"])
        sns.scatterplot(ax=axes[i, j], data=dataset_VAF, x=col, y='MutationNrPlasma')
        axes[i, j].title.set_text("Correlation p=" + str(np.round(vaf_corr[1], 4)))
        if j == 2:
            i = 1
            j = 0
        else:
            j = j + 1
    plt.subplots_adjust(bottom=0.2, right=0.8, top=0.9, wspace=0.8)
    save_file = "/Users/alexandra/PhD/FragmentationPatterns/Summary/Mutations/Plots/mutations_and_labels_scatterplot.pdf"
    plt.savefig(save_file)
    plt.clf()

    fig, axes = plt.subplots(2, 4, figsize=(18, 10))
    i = 0
    j = 0
    print("Create plots for no of mutations in plasma and  and frag sizes counts")
    for col in summary_cols_scatterplot_frag_size:
        fig.suptitle("Create plots for no of mutations in plasma and  and frag sizes counts")
        dataset_VAF = all_data[["MutationNrPlasma", col]].dropna()
        vaf_corr = stats.pearsonr(dataset_VAF[col], dataset_VAF["MutationNrPlasma"])
        sns.scatterplot(ax=axes[i, j], data=dataset_VAF, x=col, y='MutationNrPlasma')
        axes[i, j].title.set_text("Correlation p=" + str(np.round(vaf_corr[1], 4)))
        if j == 3:
            i = 1
            j = 0
        else:
            j = j + 1
    plt.subplots_adjust(bottom=0.2, right=0.8, top=0.9, wspace=0.8)
    save_file = "/Users/alexandra/PhD/FragmentationPatterns/Summary/Mutations/Plots/mutations_and_frag_counts_scatterplot.pdf"
    plt.savefig(save_file)
    plt.clf()

    fig, axes = plt.subplots(1, 3, figsize=(18, 10))
    j = 0
    print("Create plots for no of mutations in plasma and categorical labels")
    for col in summary_cols_boxplot:
        fig.suptitle("Create plots for no of mutations in plasma and categorical labels")
        dataset_VAF = all_data[["MutationNrPlasma", col]].dropna()
        # vaf_corr = stats.ttest_ind(dataset_VAF[col], dataset_VAF["MutationNrPlasma"])
        sns.boxplot(ax=axes[j], data=dataset_VAF, x=col, y='MutationNrPlasma')
        # axes[j].title.set_text("Correlation p=" + str(np.round(vaf_corr[1], 4)))
        if j == 2:
            j = 0
        j = j + 1
    plt.subplots_adjust(bottom=0.2, right=0.8, top=0.9, wspace=0.8)
    save_file = "/Users/alexandra/PhD/FragmentationPatterns/Summary/Mutations/Plots/mutations_and_labels_boxplot.pdf"
    plt.savefig(save_file)
    plt.clf()

def gather_TF_perStudy(study):
    inchorCNATF_files = glob.glob("/Users/alexandra/PhD/" + study + "/FragmentationPatterns/CNAs/ichorCNATF/*.txt")
    tumorFraction = []
    for ichor_tf_file in inchorCNATF_files:
        file_name = os.path.basename(ichor_tf_file).replace(".params.txt", "")
        with open(ichor_tf_file) as fh:
            for line in fh:
                if line.startswith("Tumor Fraction:"):
                    tumorFraction.append([file_name, str.replace(str.replace(line, "Tumor Fraction:\t", ""),
                                          "\n", "")])

    tf_per_study = pd.DataFrame(tumorFraction, columns=["SampleID", "ichorCNATF"])
    tf_per_study.to_csv("/Users/alexandra/PhD/FragmentationPatterns/Summary/CNAs/ichorCNATF_" + study + ".csv",
                        index=False, sep="\t")

def plot_histo_ichorCNATF():
    files = glob.glob("/Users/alexandra/PhD/FragmentationPatterns/Summary/CNAs/ichorCNATF_*.csv")
    for file in files:
        data = pd.read_csv(file)
        study = os.path.basename(file).replace(".csv", "").replace("ichorCNATF_", "")
        print("Run analysis for: " + study)
        histo_plot_output = "/Users/alexandra/PhD/FragmentationPatterns/Summary/CNAs/Plots/ichorTFhisto_" + study + ".pdf"
        title = "ichorCNA TF histogram for " + study
        sns.histplot(data=data, x="ichorCNATF")
        plt.xlabel("ichorCNA tumor fraction")
        plt.title(title, fontsize=8)
        plt.savefig(histo_plot_output)
        plt.clf()


def agg_no_cnas_per_sample(study):
    CNA_files = glob.glob("/Users/alexandra/PhD/" + study + "/FragmentationPatterns/CNAs/CNAsPerSample/*.csv")
    abberations_width = []

    for CNA_file in CNA_files:
        sampleID = os.path.basename(CNA_file).replace(".csv", "")
        cna_data = pd.read_csv(CNA_file, sep="\t")
        cna_data["width"] = cna_data["End"] - cna_data["Start"]
        cna_data["Abberated"] = np.where(cna_data["CNA"] > 0, 1, 0)
        agg_width_by_abb = cna_data.groupby(["Abberated"]).sum("width").reset_index()
        agg_width_by_abb["norm_width"] = np.round(agg_width_by_abb["width"]/agg_width_by_abb["width"].sum(),2)
        [abberations_width_samples] = agg_width_by_abb[agg_width_by_abb["Abberated"] == 1]\
            [["norm_width"]].values.tolist()[0]
        abberations_width.append([sampleID, abberations_width_samples])

    abberations_width_per_study = pd.DataFrame(abberations_width, columns=["SampleID", "GenomicFractionAbberated"])
    abberations_width_per_study.to_csv("/Users/alexandra/PhD/FragmentationPatterns/Summary/CNAs/genomicFractionAbberated_"
                                       + study + ".csv", index=False, sep="\t")

def plot_histo_ichorAbberations():
    files = glob.glob("/Users/alexandra/PhD/FragmentationPatterns/Summary/CNAs/genomicFractionAbberated_*.csv")
    for file in files:
        data = pd.read_csv(file, sep="\t")
        study = os.path.basename(file).replace(".csv", "").replace("genomicFractionAbberated_", "")
        print("Run analysis for: " + study)
        histo_plot_output = "/Users/alexandra/PhD/FragmentationPatterns/Summary/CNAs/Plots/" \
                            "genomicFractionAbberatedHisto_" + study + ".pdf"
        title = "ichorCNA genomic fraction aberrated histogram for " + study
        sns.histplot(data=data, x="GenomicFractionAbberated")
        plt.xlabel("genomic fraction aberrated by ichorCNA")
        plt.title(title, fontsize=8)
        plt.savefig(histo_plot_output)
        plt.clf()


def agg_tumor_vs_sample_common_abb():
    CNA_files = glob.glob("/Users/alexandra/PhD/NeoRheaStudy/FragmentationPatterns/CNAs/CNAsPerSample/*.csv")
    common_abberations_width = []
    meta_data = pd.read_csv("/Users/alexandra/PhD/NeoRheaStudy/FragmentationPatterns/Neorhea_INIVATA_results.csv",
                            sep="\t")
    for CNA_file in CNA_files:
        sampleID = os.path.basename(CNA_file).replace(".csv", "")
        meta_data_per_patient = meta_data[meta_data["SampleID"] == sampleID]
        if not meta_data_per_patient["mean_VAF"].isna().tolist()[0]:
            cna_data = pd.read_csv(CNA_file, sep="\t")
            cna_data_aberrated = cna_data[cna_data["CNA"] > 0]
            cna_data_aberrated["width"] = cna_data_aberrated["End"] - cna_data_aberrated["Start"]

            agg_width_by_abb = cna_data_aberrated.groupby(["CNA"]).sum("width").reset_index()
            agg_width_by_abb["norm_width"] = np.round(agg_width_by_abb["width"] / agg_width_by_abb["width"].sum(), 2)
            abberations_width_samples = agg_width_by_abb[agg_width_by_abb["CNA"] == 2][["norm_width"]].values.tolist()
            if len(abberations_width_samples) > 0:
                [genomic_fraction_aberrated] = abberations_width_samples[0]
            else:
                genomic_fraction_aberrated = 0
            common_abberations_width.append([sampleID, genomic_fraction_aberrated])

    abberations_width_per_study = pd.DataFrame(common_abberations_width, columns=["SampleID", "GenomicFractionAbberated"])
    abberations_width_per_study.to_csv(
        "/Users/alexandra/PhD/FragmentationPatterns/Summary/CNAs/commonGenomicFractionAbberated_NeoRheaStudy.csv",
        index=False, sep="\t")


def plot_histo_common_ichorAbberations():
    data = pd.read_csv("/Users/alexandra/PhD/FragmentationPatterns/Summary/CNAs/"
                       "commonGenomicFractionAbberated_NeoRheaStudy.csv", sep="\t")
    study = "Neorhea"
    print("Run analysis for: " + study)
    histo_plot_output = "/Users/alexandra/PhD/FragmentationPatterns/Summary/CNAs/Plots/" \
                        "commonGenomicFractionAbberatedHisto_" + study + ".pdf"
    title = "ichorCNA tumor and plasma common genomic fraction aberrated histogram for " + study
    sns.histplot(data=data, x="GenomicFractionAbberated")
    plt.xlabel("common tumor and plasma genomic fraction aberrated by ichorCNA")
    plt.title(title, fontsize=8)
    plt.savefig(histo_plot_output)
    plt.clf()


def plot_study_prop_corr_with_cnas(study):
    neorhea_data = pd.read_csv("/Users/alexandra/PhD/FragmentationPatterns/Summary/" + study + "_fragCounts.csv",
                               sep="\t")
    mutations_data = pd.read_csv("/Users/alexandra/PhD/FragmentationPatterns/Summary/CNAs/genomicFractionAbberated_"
                                 + study + ".csv", sep="\t")
    neorhea_data.rename(columns={"PatientId": "SampleID"}, inplace=True)
    all_data = pd.merge(neorhea_data, mutations_data, how="inner", on=["SampleID"])
    summary_cols_scatterplot_labels = ["ichorTF", "VAF", "TotalNrFragmentRaw", "TotalNrFragmentDELFIRaw",
                                       "TotalNrFragmentGCCorrected", "TotalNrFragmentDELFIGCCorrected"]
    summary_cols_scatterplot_frag_size = ["very_short_30_100_GC", "short_n_100_150_GC",
                                "long_n_151_220_GC", "short_din_221_300_GC", "long_din_301_360_GC",
                                "short_trin_361_500_GC", "long_trin_501_561_GC", "very_long_561_700_GC"]
    summary_cols_boxplot = ["ctDNADetected", "survivalStatus", "timepoint"]
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    i = 0
    j = 0
    print("Create plots for fraction of aberrated genomic regions in plasma and continuous labels")
    for col in summary_cols_scatterplot_labels:
        fig.suptitle(print("Create plots for fraction of aberrated genomic regions in plasma and continuous labels"))
        dataset_VAF = all_data[["GenomicFractionAbberated", col]].dropna()
        vaf_corr = stats.pearsonr(dataset_VAF[col], dataset_VAF["GenomicFractionAbberated"])
        sns.scatterplot(ax=axes[i, j], data=dataset_VAF, x=col, y='GenomicFractionAbberated')
        axes[i, j].title.set_text("Correlation p=" + str(np.round(vaf_corr[1], 4)))
        if j == 2:
            i = 1
            j = 0
        else:
            j = j + 1
    plt.subplots_adjust(bottom=0.2, right=0.8, top=0.9, wspace=0.8)
    save_file = "/Users/alexandra/PhD/FragmentationPatterns/Summary/CNAs/Plots/cnas_and_labels_scatterplot.pdf"
    plt.savefig(save_file)
    plt.clf()

    fig, axes = plt.subplots(2, 4, figsize=(18, 10))
    i = 0
    j = 0
    print("Create plots for fraction of aberrated genomic regions in plasma and  and frag sizes counts")
    for col in summary_cols_scatterplot_frag_size:
        fig.suptitle("Create plots for fraction of aberrated genomic regions in plasma and  and frag sizes counts")
        dataset_VAF = all_data[["GenomicFractionAbberated", col]].dropna()
        vaf_corr = stats.pearsonr(dataset_VAF[col], dataset_VAF["GenomicFractionAbberated"])
        sns.scatterplot(ax=axes[i, j], data=dataset_VAF, x=col, y='GenomicFractionAbberated')
        axes[i, j].title.set_text("Correlation p=" + str(np.round(vaf_corr[1], 4)))
        if j == 3:
            i = 1
            j = 0
        else:
            j = j + 1
    plt.subplots_adjust(bottom=0.2, right=0.8, top=0.9, wspace=0.8)
    save_file = "/Users/alexandra/PhD/FragmentationPatterns/Summary/CNAs/Plots/cnas_and_frag_counts_scatterplot.pdf"
    plt.savefig(save_file)
    plt.clf()

    fig, axes = plt.subplots(1, 3, figsize=(18, 10))
    j = 0
    print("Create plots for fraction of aberrated genomic regions in plasma and categorical labels")
    for col in summary_cols_boxplot:
        fig.suptitle("Create plots  for fraction of aberrated genomic regions in plasma and categorical labels")
        dataset_VAF = all_data[["GenomicFractionAbberated", col]].dropna()
        sns.boxplot(ax=axes[j], data=dataset_VAF, x=col, y='GenomicFractionAbberated')
        # axes[j].title.set_text("Correlation p=" + str(np.round(vaf_corr[1], 4)))
        if j == 2:
            j = 0
        j = j + 1
    plt.subplots_adjust(bottom=0.2, right=0.8, top=0.9, wspace=0.8)
    save_file = "/Users/alexandra/PhD/FragmentationPatterns/Summary/CNAs/Plots/cnas_and_labels_boxplot.pdf"
    plt.savefig(save_file)
    plt.clf()




def plot_response_corr_with_frag_size_2(file, output_file, cols, title):
    dataset = pd.read_csv(file, sep="\t")
    response_cols = ["Label", "ctDNADetected", "ichorTF", "survivalStatus", "timepoint"]
    if "Synergy" in title:
        dataset = dataset[~dataset["PatientId"].str.contains("Genome-IJB")]
    #check only against healthy of 50 bp
    # first_batch_samples = pd.read_csv("/Users/alexandra/PhD/healthy_sWGS/first_batch_healthy_sample_names.csv")
    # dataset = dataset[(dataset['PatientId'].isin(first_batch_samples["PatientId"].values)) | (dataset['Label'] == "Cancer")]


    dataset_cancer = dataset[dataset["Label"] == "Cancer"]
    for col in cols:
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        print("Create plot for frag sizes:" + col)
        fig.suptitle("Correlation for fragments sizes: " + col + ";" + title)
        ## Label ##
        label_test = stats.mannwhitneyu(dataset[dataset["Label"] == "Cancer"][col],
                                     dataset[dataset["Label"] == "Healthy"][col])
        sns.boxplot(ax=axes[0, 0], data=dataset, x=col, y='Label')
        axes[0, 0].title.set_text("Mann Whitney U Test p=" + str(np.round(label_test.pvalue, 4)))

        ## ctDNA detection ##
        ## some bad thing: TODO: rethink this ##
        if "Pearl" in file:
            dataset["ctDNADetected"] = np.where(dataset["VAF"] > 0, True, False)
        elif "Synergy" in file:
            dataset["ctDNADetected"] = np.where(dataset["ichorTF"] > 0.03, True, False)
        ###
        dataset_ctDNADetection = dataset[["ctDNADetected", "Label", col]].dropna()
        dataset_ctDNADetection.loc[dataset_ctDNADetection["ctDNADetected"] == True, "ctDNADetected"] = "ctDNA Detected"
        dataset_ctDNADetection.loc[dataset_ctDNADetection["ctDNADetected"] == False, "ctDNADetected"] = \
            "ctDNA not Detected"
        dataset_ctDNADetection_cancer = dataset_ctDNADetection[dataset_ctDNADetection["Label"]=="Cancer"]
        if len(dataset_ctDNADetection_cancer[
                   dataset_ctDNADetection_cancer["ctDNADetected"] == "ctDNA Detected"])>0\
                and len(dataset_ctDNADetection_cancer[
                            dataset_ctDNADetection_cancer["ctDNADetected"] == "ctDNA not Detected"])>0:

            ctDNADetection_test = stats.mannwhitneyu(
                dataset_ctDNADetection_cancer[dataset_ctDNADetection_cancer["ctDNADetected"] == "ctDNA Detected"][col],
                dataset_ctDNADetection_cancer[dataset_ctDNADetection_cancer["ctDNADetected"] == "ctDNA not Detected"][col])
            ctDNADetection_test = ctDNADetection_test.pvalue
        else:
            ctDNADetection_test = 1

        sns.boxplot(ax=axes[0, 1], data=dataset_ctDNADetection_cancer, x=col, y='ctDNADetected')
        axes[0, 1].title.set_text("Mann Whitney U Test  p=" + str(np.round(ctDNADetection_test, 4)))

        ## survival test ###
        dataset_survivalStatus = dataset[["survivalStatus", "Label", col]].dropna()
        dataset_survivalStatus.loc[dataset_survivalStatus["survivalStatus"] == 1, "survivalStatus"] = "Relapse"
        dataset_survivalStatus.loc[dataset_survivalStatus["survivalStatus"] == 0, "survivalStatus"] = "No relapse"
        dataset_survivalStatus_cancer = dataset_survivalStatus[dataset_survivalStatus["Label"] == "Cancer"]

        if len(dataset_survivalStatus_cancer[
                   dataset_survivalStatus_cancer["survivalStatus"] == "Relapse"]) > 0:
            survival_test = stats.mannwhitneyu(
                dataset_survivalStatus_cancer[dataset_survivalStatus_cancer["survivalStatus"] == "Relapse"][col],
                dataset_survivalStatus_cancer[dataset_survivalStatus_cancer["survivalStatus"] == "No relapse"][col])
            survival_test = survival_test.pvalue
        else:
            survival_test = 1
        sns.boxplot(ax=axes[0, 2], data=dataset_survivalStatus_cancer, x=col, y='survivalStatus')
        axes[0, 2].title.set_text("Mann Whitney U Test  p=" + str(np.round(survival_test, 4)))

        #### ichorTF ##
        dataset["ichorTF_threshold"] = np.where(dataset["ichorTF"] >= 0.03, "TF>=0.03", "TF<0.03")
        dataset_ichor = dataset[["ichorTF_threshold", "Label", col]].dropna()
        dataset_ichor_cancer = dataset_ichor[dataset_ichor["Label"] == "Cancer"]
        if len(dataset_ichor_cancer[
                   dataset_ichor_cancer["ichorTF_threshold"] == "TF>=0.03"])>0\
                and len(dataset_ichor_cancer[
                            dataset_ichor_cancer["ichorTF_threshold"] == "TF<0.03"])>0:

            ichor_test = stats.mannwhitneyu(
                dataset_ichor_cancer[dataset_ichor_cancer["ichorTF_threshold"] == "TF>=0.03"][col],
                dataset_ichor_cancer[dataset_ichor_cancer["ichorTF_threshold"] == "TF<0.03"][col])
            ichor_test = ichor_test.pvalue
        else:
            ichor_test = 1

        sns.boxplot(ax=axes[1, 0], data=dataset_ichor_cancer, x=col, y='ichorTF_threshold')
        axes[1, 0].title.set_text("Mann Whitney U Test  p=" + str(np.round(ichor_test, 4)))

        ### VAF ##

        # dataset_VAF = dataset[["Label", "VAF", col]].dropna()
        # dataset_VAF_cancer = dataset_cancer[["Label", "VAF", col]].dropna()
        # if dataset_VAF_cancer.empty:
        #     ichortf_corr = stats.pearsonr(dataset_cancer[col], dataset_cancer["ichorTF"])
        #     sns.scatterplot(ax=axes[1, 0], data=dataset, x=col, y='ichorTF', hue="Label")
        #     axes[1, 1].title.set_text("Correlation p=" + str(np.round(ichortf_corr[1], 4)))
        # else:
        #     vaf_corr = stats.pearsonr(dataset_VAF_cancer[col], dataset_VAF_cancer["VAF"])
        #     sns.scatterplot(ax=axes[1, 1], data=dataset_VAF, x=col, y='VAF', hue="Label")
        #     axes[1, 1].title.set_text("Correlation p=" + str(np.round(vaf_corr[1], 4)))

        ## timepoint ###
        timepoint_test = stats.kruskal(dataset_cancer[dataset_cancer["timepoint"] == "BL"][col],
                                       dataset_cancer[dataset_cancer["timepoint"] == "C1"][col],
                                       dataset_cancer[dataset_cancer["timepoint"] == "Surgery"][col])
        sns.boxplot(ax=axes[1, 1], data=dataset_cancer, x=col, y='timepoint')
        axes[1, 1].title.set_text("kruskal wallis test p=" + str(np.round(timepoint_test.pvalue, 4)))
        # for ax in axes:
        #     for ax_line in ax:
        #         # ax_line.text(0.85, 0.85, ctDNADetection_test[1], fontsize=9)

        plt.subplots_adjust(bottom=0.2, right=0.8, top=0.9, wspace=0.8)
        save_file = output_file + "_" + col + ".pdf"
        plt.savefig(save_file)
        plt.clf()