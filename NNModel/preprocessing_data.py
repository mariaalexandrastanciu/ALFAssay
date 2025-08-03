# Created by alexandra at 18/12/2023
import numpy as np
import pandas as pd
import CONSTANTS as c
import pyranges as pr
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import kstest
import glob
import os.path
from sklearn.preprocessing import PowerTransformer
import utils as u

def get_feature_columns():
    feature_columns = ["noReads", "noReadsGCParagonFrag", "ultraShortFrag", "ultraShortFragGCPara",
                       "shortFrag", "shortFragGCPara", "longFrag", "longFragGCPara", "ultraLongFrag",
                       "ultraLongFragGCPara", "shortOverLong", "shortOverLongGCP    ara"]

    return feature_columns

def create_data(file, metadata_file, output):
    data = pd.read_csv(file, sep="\t")
    metadata = pd.read_csv(metadata_file, sep="\t")

    melt_short_reads_columns = ["short_arm_" + str(x) for x in c.armlevels]
    # melt_short_reads_columns.extend(["PatientId"])

    melt_no_reads_columns = ["NoReads_arm_" + str(x) for x in c.armlevels]
    # melt_no_reads_columns.extend(["PatientId"])

    # data_short_reads = data[melt_short_reads_columns]
    data_short_reads_melted = pd.melt(data, id_vars=[("PatientId")], value_vars=[*melt_short_reads_columns],
                              var_name='region', value_name='short_reads')
    data_short_reads_melted["region"] = data_short_reads_melted["region"].str.replace("short_arm_", "")
    # data_no_reads = data[melt_no_reads_columns]
    data_no_reads_melted = pd.melt(data, id_vars=[("PatientId")],
                                      value_vars=[*melt_no_reads_columns],
                                      var_name='region', value_name='no_reads')
    data_no_reads_melted["region"] = data_no_reads_melted["region"].str.replace("NoReads_arm_", "")

    melted_data = pd.merge(data_short_reads_melted, data_no_reads_melted, on=["PatientId", "region"])
    melted_data = pd.merge(metadata, melted_data, on="PatientId")

    melted_data.to_csv(output, sep="\t", index=False)


#
# metadata_file = "/Users/alexandra/PhD/FragmentationPatterns/Data/MetaData/AllStudiesMetaData.csv"
# file = "/Users/alexandra/PhD/FragmentationPatterns/Data/version1/window_fragsize_short_noReadsByArmNoScale.csv"
# output = "/Users/alexandra/PhD/FragmentationPatterns/Data/presentation_input2/window_fragsize_short_noReads_byArmNoScale.csv"
# create_data(file, metadata_file, output)
#
# def get_data_for_wrapper(file, test_studies, validation_studies, task ):
#     data = pd.read_csv(file, sep="\t")
#     data = data.dropna(subset=[c.task_dict.get(task)])
#     # data = data[adta["timepoint"]=="BL"]
#     data["Label"] = np.where(data["Label"] == "Cancer", 1, 0)
#     data["ctDNADetected"] = np.where(data["ctDNADetected"] == True, 1, data["ctDNADetected"])
#     data["ctDNADetected"] = np.where(data["ctDNADetected"] == False, 0, data["ctDNADetected"])
#     data["ctDNADetected"] = np.where(data["ctDNADetected"].isin([0, 1]), data["ctDNADetected"], -1)
#
#     data["VAFg0p001"] = np.where(data["VAFg0p001"] == True, 1, data["VAFg0p001"])
#     data["VAFg0p001"] = np.where(data["VAFg0p001"] == False, 0, data["VAFg0p001"])
#     data["VAFg0p001"] = np.where(data["VAFg0p001"].isin([0, 1]), data["VAFg0p001"], -1)
#
#     test_data = data[data["study"].isin(test_studies)]
#     validation_data = data[data["study"].isin(validation_studies)]
#     # features_columns = ["short_reads", "no_reads"]
#
#     features_columns = ["noReadsGCParagonFrag", "shortOverLongGCPara"]
#
#     # features_columns = ["shortOverLongGCPara"]
#     # features_columns = ["noReadsGCParagonFrag",  "ultraShortFragGCPara", "shortFragGCPara", "longFragGCPara",
#     #                     "ultraLongFragGCPara", "shortOverLongGCPara"]
#     #
#     # features_columns = ["noReads", "noReadsGCParagonFrag", "ultraShortFrag", "ultraShortFragGCPara",
#     #                     "shortFrag", "shortFragGCPara", "longFrag", "longFragGCPara", "ultraLongFrag",
#     #                     "ultraLongFragGCPara", "shortOverLong", "shortOverLongGCPara"]
#
#     # labels_columns = ["PatientId", "Label", "ctDNADetected", "VAFg0p001"]
#
#     test_labels = test_data[c.labels_columns].drop_duplicates().values
#     fragmentSizeSampleTest = []
#     for i, sample in enumerate(test_labels):
#         features_per_sample = test_data[test_data["PatientId"]==sample[0]][features_columns]
#         test_labels_per_sample=test_labels[i]
#         fragmentSizeSampleTest.append([features_per_sample, test_labels_per_sample])
#
#     validation_labels = validation_data[c.labels_columns].drop_duplicates().values
#     fragmentSizeSampleValidation = []
#     for i, sample in enumerate(validation_labels):
#         features_per_sample = validation_data[validation_data["PatientId"]==sample[0]][features_columns]
#         validation_labels_per_sample=validation_labels[i]
#         fragmentSizeSampleValidation.append([features_per_sample, validation_labels_per_sample])
#     fragmentSizeFeatureTest = [item[0] for item in fragmentSizeSampleTest]
#     fragmentSizeFeatureValidation = [item[0] for item in fragmentSizeSampleValidation]
#
#     return fragmentSizeSampleTest, fragmentSizeSampleValidation, test_labels, validation_labels, fragmentSizeFeatureTest, \
#            fragmentSizeFeatureValidation


# file = "/Users/alexandra/PhD/FragmentationPatterns/Data/NN_5MB/NN_5MB_scaled_2.csv"
# test_studies = ["PearlStudy"]
# validation_studies = ["NeoRheaStudy", "healthy_sWGS"]
# task = 2
# fragmentSizeSampleTest, fragmentSizeSampleValidation,\
#     test_labels, validation_labels = get_data_for_wrapper(file, test_studies, validation_studies,
#                                                                             task)
# print("ok")
# #

def add_DELFI_metadata():
    meta_data = pd.read_csv("/Users/alexandra/PhD/FragmentationPatterns/Data/MetaData/AllStudiesMetaData.csv", sep="\t")
    meta_data_delfi = pd.read_csv("/Users/alexandra/PhD/FragmentationPatterns/Data/MetaData/DELFIStudyMetaData.csv", sep="\t")
    meta_data_delfi["study"] = "DELFI"

    result = pd.concat([meta_data, meta_data_delfi])

    result.to_csv("/Users/alexandra/PhD/FragmentationPatterns/Data/MetaData/AllStudiesMetaDataWDELFI.csv", sep="\t",
                  index=False)
# add_DELFI_metadata()

def norm_per_sample(sample_data):
    data_scaled = sample_data[["PatientId", "bin", "Chromosome", "Start", "End", "arm", "gcContent",
                               "ultraShortFrag", "ultraShortFragGCPara",
                                                  "shortFrag", "shortFragGCPara", "longFrag", "longFragGCPara", "ultraLongFrag",
                                                  "ultraLongFragGCPara",  "shortOverLong", "shortOverLongGCPara"
                               ]]
    # feature_columns = ["noReads", "noReadsGCParagonFrag", "ultraShortFrag", "ultraShortFragGCPara",
    #                    "shortFrag", "shortFragGCPara", "longFrag", "longFragGCPara", "ultraLongFrag",
    #                    "ultraLongFragGCPara",]

    feature_columns = ["noReads", "noReadsGCParagonFrag"]
    for col in feature_columns:
        data_scaled[col] = u.basic_scale(sample_data[col])

    return data_scaled


def norm_per_sample2(sample_data):
    data_scaled = sample_data[["PatientId", "bin", "Chromosome", "Start", "End", "arm", "gcContent",
                               "ultraShortFrag", "ultraShortFragGCPara",
                                                  "shortFrag", "shortFragGCPara", "longFrag", "longFragGCPara", "ultraLongFrag",
                                                  "ultraLongFragGCPara",  "shortOverLong", "shortOverLongGCPara"
                               ]]
    # feature_columns = ["noReads", "noReadsGCParagonFrag", ]
    feature_columns = ["noReads", "noReadsGCParagonFrag"]
    for col in feature_columns:
        data_scaled[col] = u.norm_alternative(sample_data[col])

    return data_scaled

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

# studies = ["healthy_sWGS", "NeoRheaStudy", "PearlStudy", "SynergyStudy"]
# # studies = ["healthy_sWGS"]
# # studies = ["SynergyStudy"]
# root = "/Users/alexandra/PhD/"
# first_file = True
# write_mode = "w"
# analysis_input = "NN_5MB" #"WindowFragmentSizes_30_700_mq40_wz100k" #"WindowFragmentSizes_hg38_30_700_mq60" #"WindowFragmentSizes_30_700_mq40_wz100k"
# analysis_output = "NN_5MB"
# # meta_data = pd.read_csv("/Users/alexandra/PhD/FragmentationPatterns/Data/MetaData/AllStudiesMetaDataWDELFI.csv", sep="\t")
# for i, study in enumerate(studies):
#     if study == "healthy_sWGS":
#         meta = "Healthy"
#     else:
#         meta = study
#     if i > 0:
#         first_file = False
#         write_mode = "a"
#     meta_data = pd.read_csv(root + "FragmentationPatterns/Data/MetaData/" + meta + "MetaData.csv", sep="\t")
#     # meta_data = meta_data[meta_data["Label"] == "Healthy_Synergy"]
#
#     input_files = sorted(glob.glob(root + study +
#                             "/FragmentationPatterns/FragmentSizes/" + analysis_input + "/*.bed"))
#     output_folder = root + "FragmentationPatterns/Data/" + analysis_output + "/NN_5MB.csv"
#     output_folder_scaled = root + "FragmentationPatterns/Data/" + analysis_output + "/NN_5MB_scaled.csv"
#     output_folder_scaled_2 = root + "FragmentationPatterns/Data/" + analysis_output + "/NN_5MB_scaled_2.csv"
#     concat_samples(input_files, output_folder, meta_data, study, first_file, write_mode, output_folder_scaled,
#                    output_folder_scaled_2)
# #
def concat_input_at_read_level(sample_data_file, meta_data, study, output_file, write_mode, first_file):
    sample_data = pd.read_csv(sample_data_file, sep="\t")
    features = sample_data[["fragment_size", "GC", "gcParagonWeight", "end4BP"]]
    sample_id = u.get_patient_from_file_name(sample_data_file, "_fragment_size_expanded_s.bed")
    meta_data_per_patient = meta_data[meta_data["PatientId"] == sample_id]
    # data_per_patient = pd.concat([features, meta_data_per_patient])
    features["PatientId"] = sample_id
    data_per_patient = pd.merge(features, meta_data_per_patient, on="PatientId")
    data_per_patient["study"] = study
    data_per_patient.to_csv(output_file, sep="\t", mode=write_mode, index=False, header=first_file)
    return data_per_patient

# root = "/Users/alexandra/PhD/"
# studies = ["healthy_sWGS", "PearlStudy"]
# output_file = "/Users/alexandra/PhD/FragmentationPatterns/Data/NN_5MB/NN_5MB_reads.csv"
# first_file = True
# write_mode = "w"
#
# for i, study in enumerate(studies):
#     path = "/Users/alexandra/PhD/" + study + "/FragmentationPatterns/FragmentSizes/expanded/"
#     files = glob.glob(path + "*expanded_s.bed")
#     if study == "healthy_sWGS":
#         meta = "healthy"
#         meta_data = pd.read_csv(root + "FragmentationPatterns/Data/MetaData/" + meta + "MetaData.csv", sep="\t")
#     else:
#         meta_data = pd.read_csv(root + "FragmentationPatterns/Data/MetaData/" + study + "MetaData.csv", sep="\t")
#
#     if i > 0:
#         first_file = False
#         write_mode = "a"
#     for file in files:
#         concat_input_at_read_level(file, meta_data, study, output_file, write_mode, first_file)

def frag_distribution(sample_data_file, output_file):
    sample_data = pd.read_csv(sample_data_file, sep="\t")
    fragment_size_distribution = sample_data.groupby(['fragment_size'])['gcParagonWeight'].sum()
    fragment_size_distribution.to_csv(output_file, sep="\t")

# sample_data_file = "/Users/alexandra/PhD/healthy_sWGS/FragmentationPatterns/FragmentSizes/expanded/Genome-IJB-HP-9-xx_S39_fragment_size_expanded_sr.bed"
# output_file = "/Users/alexandra/PhD/healthy_sWGS/FragmentationPatterns/FragmentSizes/frag_distrib/Genome-IJB-HP-9-xx_S39_fragment_ditrib_sr.csv"
# frag_distribution(sample_data_file, output_file)


def visualize_density(sample_data_file_1, sample_data_file_2, plot_output):
    sample_data_1 = pd.read_csv(sample_data_file_1, sep="\t")
    sample_data_2 = pd.read_csv(sample_data_file_2, sep="\t")
    # sample_data_1["gcParagonWeight"] = (sample_data_1["gcParagonWeight"] - sample_data_1["gcParagonWeight"].mean())/sample_data_1["gcParagonWeight"].std()
    # sample_data_2["gcParagonWeight"] = (sample_data_2["gcParagonWeight"] - sample_data_2["gcParagonWeight"].mean()) / \
    #                                    sample_data_2["gcParagonWeight"].std()
    # sample_data_1["gcParagonWeight"] = (1e6*sample_data_1["gcParagonWeight"]/sample_data_1["gcParagonWeight"].sum())
    # sample_data_2["gcParagonWeight"] = (1e6 * sample_data_2["gcParagonWeight"] / sample_data_2["gcParagonWeight"].sum())
    sample_data_1["Label"] = "Healthy"
    sample_data_2["Label"] = "Pearl"
    all_data = pd.concat([sample_data_1, sample_data_2])

    sns.lineplot(x="fragment_size", y="gcParagonWeight", hue="Label", data=all_data)
    plt.xlabel("Fragment Sizes Distribution")
    plt.title("Healthy vs Pearl all regions", fontsize=8)
    plt.savefig(plot_output)
    plt.clf()


# sample_data_file_1 = "/Users/alexandra/PhD/healthy_sWGS/FragmentationPatterns/FragmentSizes/frag_distrib/Genome-IJB-HP-9-xx_S39_fragment_ditrib.csv"
# sample_data_file_2 = "/Users/alexandra/PhD/PearlStudy/FragmentationPatterns/FragmentSizes/frag_distrib/NIPT-PRL-049-25-BL-Px_S61_fragment_ditrib.csv"
# plot_output = "/Users/alexandra/PhD/healthy_sWGS/FragmentationPatterns/FragmentSizes/plots/healthy_vs_pearl.pdf"
# visualize_density(sample_data_file_1, sample_data_file_2, plot_output)


def plot_frag_size_on_motif(sample1_file, sample2_file,plot_output, title, motif):
    sample_data_1 = pd.read_csv(sample1_file, sep="\t")
    sample_data_2 = pd.read_csv(sample2_file, sep="\t")

    if motif:
        # sample_data_1 = sample_data_1[sample_data_1["end4BP"]=="CCCA"]
        sample_data_2 = sample_data_2[sample_data_2["end4BP"]=="CCCA"]

    fragment_size_distribution_1 = sample_data_1.groupby(['fragment_size'])['gcParagonWeight'].sum().reset_index()
    fragment_size_distribution_2 = sample_data_2.groupby(['fragment_size'])['gcParagonWeight'].sum().reset_index()

    fragment_size_distribution_1["gcParagonWeight"] = u.basic_scale(fragment_size_distribution_1["gcParagonWeight"])
    fragment_size_distribution_2["gcParagonWeight"] = u.basic_scale(fragment_size_distribution_2["gcParagonWeight"])

    if motif:
        fragment_size_distribution_1 = fragment_size_distribution_1[fragment_size_distribution_1["end4BP"]=="CCCA"]
        fragment_size_distribution_2 = fragment_size_distribution_2[fragment_size_distribution_2["end4BP"]=="CCCA"]

    fragment_size_distribution_1["Label"] = "Healthy_specific"
    fragment_size_distribution_2["Label"] = "Healthy_specific_CCCA"
    all_data = pd.concat([fragment_size_distribution_1, fragment_size_distribution_2])

    test_result = kstest(fragment_size_distribution_1["gcParagonWeight"].values,
                         fragment_size_distribution_2["gcParagonWeight"].values, alternative='less')

    print(test_result)
    sns.lineplot(x="fragment_size", y="gcParagonWeight", hue="Label", data=all_data)
    plt.xlabel("Fragment Sizes Distribution")
    title = title + " KV test " + str(np.round(test_result.statistic)) + ", p:" + str(np.round(test_result.pvalue))

    plt.title(title, fontsize=6)
    plt.savefig(plot_output)
    plt.clf()


regions = ["all", "specific"]
motif=True
root ="/Users/alexandra/PhD/FragmentationPatterns/Summary/Plots/Enrichment/"

# if motif:
#     for region in regions:
#         if region == "all":
#             sample_data_file_1 = "/Users/alexandra/PhD/healthy_sWGS/FragmentationPatterns/FragmentSizes/expanded/Genome-IJB-HP-9-xx_S39_fragment_size_expanded.bed"
#             sample_data_file_2 = "/Users/alexandra/PhD/PearlStudy/FragmentationPatterns/FragmentSizes/expanded/NIPT-PRL-049-25-BL-Px_S61_fragment_size_expanded.bed"
#             plot_output = root + "pearl_healthy_all_regions_CCCA_scaled.pdf"
#             title = "Fragment Sizes Distribution cancer vs healthy all regions, CCCA motif"
#             plot_frag_size_on_motif(sample_data_file_1, sample_data_file_2, plot_output, title, motif)
#         else:
#             sample_data_file_1 = "/Users/alexandra/PhD/healthy_sWGS/FragmentationPatterns/FragmentSizes/expanded/Genome-IJB-HP-9-xx_S39_fragment_size_expanded_sr.bed"
#             sample_data_file_2 = "/Users/alexandra/PhD/PearlStudy/FragmentationPatterns/FragmentSizes/expanded/NIPT-PRL-049-25-BL-Px_S61_fragment_size_expanded_sr.bed"
#             plot_output = root + "pearl_healthy_specific_regions_CCCA_scaled.pdf"
#             title = "Fragment Sizes Distribution cancer vs healthy specific regions, CCCA motif"
#             plot_frag_size_on_motif(sample_data_file_1, sample_data_file_2, plot_output, title, motif)
# else:
#     for region in regions:
#         if region == "all":
#             sample_data_file_1 = "/Users/alexandra/PhD/healthy_sWGS/FragmentationPatterns/FragmentSizes/expanded/Genome-IJB-HP-9-xx_S39_fragment_size_expanded.bed"
#             sample_data_file_2 = "/Users/alexandra/PhD/PearlStudy/FragmentationPatterns/FragmentSizes/expanded/NIPT-PRL-049-25-BL-Px_S61_fragment_size_expanded.bed"
#             plot_output =  root + "pearl_healthy_all_regions_scaled.pdf"
#             title = "Fragment Sizes Distribution cancer vs healthy all regions"
#             plot_frag_size_on_motif(sample_data_file_1, sample_data_file_2, plot_output, title, motif)
#         else:
#             sample_data_file_1 = "/Users/alexandra/PhD/healthy_sWGS/FragmentationPatterns/FragmentSizes/expanded/Genome-IJB-HP-9-xx_S39_fragment_size_expanded_sr.bed"
#             sample_data_file_2 = "/Users/alexandra/PhD/PearlStudy/FragmentationPatterns/FragmentSizes/expanded/NIPT-PRL-049-25-BL-Px_S61_fragment_size_expanded_sr.bed"
#             plot_output =  root + "pearl_healthy_specific_regions_scaled.pdf"
#             title = "Fragment Sizes Distribution cancer vs healthy specific regions"
#             plot_frag_size_on_motif(sample_data_file_1, sample_data_file_2, plot_output, title, motif)


# sample_data_file_1 = "/Users/alexandra/PhD/healthy_sWGS/FragmentationPatterns/FragmentSizes/expanded/Genome-IJB-HP-9-xx_S39_fragment_size_expanded_sr.bed"
# sample_data_file_2 = "/Users/alexandra/PhD/healthy_sWGS/FragmentationPatterns/FragmentSizes/expanded/Genome-IJB-HP-9-xx_S39_fragment_size_expanded_sr.bed"
# plot_output = root + "healthy_specific_regions_CCCA_scaled.pdf"
# title = "Fragment Sizes Distribution healthy specific regions, CCCA motif"
# plot_frag_size_on_motif(sample_data_file_1, sample_data_file_2, plot_output, title, True)

# sample_data_file_1 = "/Users/alexandra/PhD/healthy_sWGS/FragmentationPatterns/FragmentSizes/expanded/Genome-IJB-HP-9-xx_S39_fragment_size_expanded.bed"
# sample_data_file_2 = "/Users/alexandra/PhD/healthy_sWGS/FragmentationPatterns/FragmentSizes/expanded/Genome-IJB-HP-9-xx_S39_fragment_size_expanded.bed"
# plot_output = root + "healthy_specific_regions_scaled.pdf"
# title = "Fragment Sizes Distribution, healthy all regions, CCCA motif"
# plot_frag_size_on_motif(sample_data_file_1, sample_data_file_2, plot_output, title, True)

# sample_data_file_1 = "/Users/alexandra/PhD/healthy_sWGS/FragmentationPatterns/FragmentSizes/expanded/Genome-IJB-HP-9-xx_S39_fragment_size_expanded.bed"
# sample_data_file_2 = "/Users/alexandra/PhD/healthy_sWGS/FragmentationPatterns/FragmentSizes/expanded/Genome-IJB-HP-9-xx_S39_fragment_size_expanded_sr.bed"
# plot_output = root + "healthy_all_specufic_regions_scaled.pdf"
# title = "Fragment Sizes Distribution healthy all vs specific regions"
# plot_frag_size_on_motif(sample_data_file_1, sample_data_file_2, plot_output, title, False)


# sample_data_file_1 = "/Users/alexandra/PhD/healthy_sWGS/FragmentationPatterns/FragmentSizes/expanded/Genome-IJB-HP-9-xx_S39_fragment_size_expanded.bed"
# sample_data_file_2 = "/Users/alexandra/PhD/healthy_sWGS/FragmentationPatterns/FragmentSizes/expanded/Genome-IJB-HP-9-xx_S39_fragment_size_expanded_sr.bed"
# plot_output = root + "healthy_all_specufic_regions_CCCA_scaled.pdf"
# title = "Fragment Sizes Distribution healthy all vs specific regions, CCCA motif"
# plot_frag_size_on_motif(sample_data_file_1, sample_data_file_2, plot_output, title, True)


def create_amplified_reads(sample_data_file, amp_regions_file, output_file):
    sample_data = pr.PyRanges(pd.read_csv(sample_data_file, sep="\t"))
    specific_regions_pr = pr.PyRanges(pd.read_csv(amp_regions_file, sep="\t"))
    result = sample_data.overlap(specific_regions_pr).as_df()
    result.to_csv(output_file, sep="\t")

# sample_file_data="/Users/alexandra/PhD/healthy_sWGS/FragmentationPatterns/FragmentSizes/expanded/Genome-IJB-HP-9-xx_S39_fragment_size_expanded.bed"
# output_file="/Users/alexandra/PhD/healthy_sWGS/FragmentationPatterns/FragmentSizes/expanded/Genome-IJB-HP-9-xx_S39_fragment_size_expanded_sr_s.bed"
# amp_regions_file = "/Users/alexandra/PhD/PyCharmProjects/ALFAssay/refdata/amp_segments.bed"
# create_amplified_reads(sample_file_data, amp_regions_file, output_file)
# path = "/Users/alexandra/PhD/PearlStudy/FragmentationPatterns/FragmentSizes/expanded/"
# files = glob.glob(path + "*NIPT-PRL-060-17-PD-Px_S129_fragment_size_expanded.bed")
# amp_regions_file = "/Users/alexandra/PhD/PyCharmProjects/ALFAssay/refdata/amp_segments.bed"
# for file in files:
#     output_file = os.path.basename(file).replace(".bed", "_s.bed")
#     create_amplified_reads(file, amp_regions_file, path+output_file)



def create_read_input_nn(file, test_studies, validation_studies, task ):
    data = pd.read_csv(file, sep="\t")
    data = data.dropna(subset=[c.task_dict.get(task)])

    data["Label"] = np.where(data["Label"] == "Cancer", 1, 0)
    data["ctDNADetected"] = np.where(data["ctDNADetected"] == True, 1, data["ctDNADetected"])
    data["ctDNADetected"] = np.where(data["ctDNADetected"] == False, 0, data["ctDNADetected"])
    data["ctDNADetected"] = np.where(data["ctDNADetected"].isin([0, 1]), data["ctDNADetected"], -1)

    data["VAFg0p001"] = np.where(data["VAFg0p001"] == True, 1, data["VAFg0p001"])
    data["VAFg0p001"] = np.where(data["VAFg0p001"] == False, 0, data["VAFg0p001"])
    data["VAFg0p001"] = np.where(data["VAFg0p001"].isin([0, 1]), data["VAFg0p001"], -1)

    test_data = data[data["study"].isin(test_studies)]
    validation_data = data[data["study"].isin(validation_studies)]

    features_columns = c.read_features

    labels_columns = c.labels_columns

    test_labels_subpsample = test_data.groupby("ctDNADetected").sample(n=10000, random_state=1).dropna()
    fragmentSizeFeatureTest = test_labels_subpsample[features_columns]
    test_labels = test_labels_subpsample[labels_columns]
    
    return test_labels_subpsample, fragmentSizeFeatureTest, test_labels


def gc_content_check():
    data = pd.read_csv("/Users/alexandra/PhD/FragmentationPatterns/Data/NN_5MB/NN_5MB_scaled.csv", sep="\t")
    data = data[data["PatientId"]=="Genome-IJB-HP-1-xx_S73"]
    pair_plot_data=data[["gcContent", "noReadsGCParagonFrag", "noReads", "shortOverLongGCPara", "shortOverLong"]]
    sns.pairplot(pair_plot_data)
    plt.title("Correlation plot 5MB", fontsize=8)
    plt.savefig("/Users/alexandra/PhD/GCCorrection/gccorrectioncheck/correlation_plot_5MB_patient.pdf")
    plt.clf()

    # data_read = pd.read_csv("/Users/alexandra/PhD/FragmentationPatterns/Data/NN_5MB/NN_5MB_reads.csv", sep="\t")

# gc_content_check()

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

fields = ["noReadsGCParagonFrag", "shortFragGCPara", "longFragGCPara", "shortOverLongGCPara", "ultraLongFragGCPara"]
root = "/Users/alexandra/PhD/FragmentationPatterns/Data/NN_5MB/"
input_file = root + "NN_5MB.csv"

output_file = root + "NN_5MB_pivot.csv"
# pivot_data(input_file, output_file, fields)

def normalize_data(file, norm_filds, output):
    data = pd.read_csv(file, sep="\t")
    excl_data = data[~data["Feature"].isin(norm_filds)]
    data = data[data["Feature"].isin(norm_filds)]

    col = "bin"
    samples = data.loc[:, data.columns != col]
    samples_nrom = samples.groupby("Feature").apply(u.norm_up_one).reset_index()
    samples_nrom["level_1"] = data["bin"]
    samples_nrom = samples_nrom.rename(columns={"level_1": "bin"})
    samples_nrom_all = pd.concat([samples_nrom, excl_data], axis=0)
    samples_nrom_all.to_csv(output,
                        sep="\t",
                        index=False)
    return samples_nrom

file = "/Users/alexandra/PhD/FragmentationPatterns/Data/NN_5MB/NN_5MB_pivot.csv"
norm_filds = ["noReadsGCParagonFrag", "shortFragGCPara", "longFragGCPara"]
output = "/Users/alexandra/PhD/FragmentationPatterns/Data/NN_5MB/NN_5MB_pivot_scaled.csv"
# normalize_data(file, norm_filds, output)

def create_data_wrapper(all_patients, data_feat, training_columns, validation_columns, metadata, labels_cols,features ):
    data_feat_cols = ["bin"] + ["Feature"] + all_patients
    data_feat = data_feat[data_feat_cols]
    all_columns = data_feat.columns
    eps = 0.001

    wrapperDataTraining = []
    labels_training = []
    features_training = []
    for patient in training_columns:
        patient_data = data_feat[["Feature", patient]]
        metadata_patient = metadata[metadata["PatientId"] == patient]
        features_p = []
        for feat in features:
            # ["noReadsGCParagonFrag", "shortFragGCPara"]
            patient_data_feat_reads = patient_data[patient_data["Feature"] == "noReadsGCParagonFrag"][patient].apply(
                pd.to_numeric).values + eps
            patient_data_feat_short_reads = patient_data[patient_data["Feature"] == "shortFragGCPara"][patient].apply(
                pd.to_numeric).values
            patient_data_feat_ultralong_reads = patient_data[patient_data["Feature"] == "ultraLongFragGCPara"][
                patient].apply(
                pd.to_numeric).values
            features_p.append([patient_data_feat_short_reads / patient_data_feat_reads,
                               patient_data_feat_ultralong_reads / patient_data_feat_reads
                               ])
        labels = metadata_patient[labels_cols].values
        wrapperDataTraining.append([np.asarray(features_p).T, labels])
        labels_training.append(labels)
        features_training.append(np.asarray(features_p).T)

    wrapperDataValidation = []
    labels_validation = []
    features_validation = []
    for patient in validation_columns:
        patient_data = data_feat[["Feature", patient]]
        metadata_patient = metadata[metadata["PatientId"] == patient]
        features_p = []
        for feat in features:
            # patient_data_feat = patient_data[patient_data["Feature"] == feat][patient].values
            patient_data_feat_reads = patient_data[patient_data["Feature"] == "noReadsGCParagonFrag"][patient].apply(
                pd.to_numeric).values + eps
            patient_data_feat_short_reads = patient_data[patient_data["Feature"] == "shortFragGCPara"][patient].apply(
                pd.to_numeric).values
            patient_data_feat_ultralong_reads = patient_data[patient_data["Feature"] == "ultraLongFragGCPara"][
                patient].apply(
                pd.to_numeric).values
            features_p.append([patient_data_feat_short_reads / patient_data_feat_reads,
                               patient_data_feat_ultralong_reads / patient_data_feat_reads
                               ])
        labels = metadata_patient[labels_cols].values
        wrapperDataValidation.append([np.asarray(features_p).T, labels])
        labels_validation.append(labels)
        features_validation.append(np.asarray(features_p).T)

    labels_training = np.asarray(labels_training).squeeze()
    features_training = np.asarray(features_training).squeeze()
    labels_validation = np.asarray(labels_validation).squeeze()
    features_validation = np.asarray(features_validation).squeeze()

    return wrapperDataTraining, wrapperDataValidation, labels_training, features_training, labels_validation, features_validation



def get_data_wrapper(file, training_studies, validation_studies, features, task):
    metadata = pd.read_csv("/Users/alexandra/PhD/FragmentationPatterns/Data/MetaData/AllMetaData.csv", sep="\t")
    labels_cols = ["PatientId", "Label", "ctDNADetected", "VAF", "ichorTF"]
    data = pd.read_csv(file, sep="\t")
    data = data.fillna(0)

    # data_feat = data[data["Feature"].isin(features)]
    data_feat = data
    metadata = metadata[metadata[task].notna()]
    metadata["ctDNADetected"] = metadata["ctDNADetected"].apply(float)
    metadata["ichorTF"] = metadata["ichorTF"].apply(float)
    metadata["ichorTF"] = np.where(metadata["ichorTF"]<0.03, 0 ,metadata["ichorTF"]*100)
    # metadata["ichorTF"] = metadata["ichorTF"] * 100
    metadata["Label"] = metadata["Label"].apply(float)
    metadata["VAF"] = metadata["VAF"].apply(float)

    ## transform nn output
    ichorTF = metadata["ichorTF"].values
    pt = PowerTransformer(method='yeo-johnson')
    ichorTF = ichorTF.reshape(-1, 1)
    pt.fit(ichorTF)
    transformed_data = pt.transform(ichorTF)
    metadata["ichorTFND"] = transformed_data

    # temporary fix see why this sample is missing
    missing = ["Genome-IJB-HP-40-xx_S28", "NIPT-PIJB-HP-73-Normal-15052023_S143"] #
    meta_v2 = metadata[~metadata["PatientId"].isin(missing)]
    # patients = list(metadata[metadata["PatientId"] != "NIPT-PRL-062-25-D14-Px_S83"]["PatientId"].values)
    patients = list(metadata["PatientId"].values)
    # patients = list(meta_v2["PatientId"].values)

    data_feat_cols = ["bin"] + ["Feature"] + patients
    data_feat = data_feat[data_feat_cols]
    all_columns = data_feat.columns
    eps = 0.001

    training_columns = list(metadata[metadata["Process"]=="Training"]["PatientId"].values)
    # for training_study in training_studies:
    #     training_columns = training_columns + list(filter(lambda x: training_study in x, all_columns))

    validation_columns = list(metadata[metadata["Process"]=="Validation"]["PatientId"].values)
    # for validation_study in validation_studies:
    #     validation_columns = validation_columns + list(filter(lambda x: validation_study in x, all_columns))

    wrapperDataTraining = []
    labels_training = []
    features_training = []
    for patient in training_columns:
        patient_data = data_feat[["Feature", patient]]
        metadata_patient = metadata[metadata["PatientId"]==patient]
        features_p = []
        for feat in features:
            # ["noReadsGCParagonFrag", "shortFragGCPara"]
            patient_data_feat_reads = patient_data[patient_data["Feature"]=="noReadsGCParagonFrag"][patient].apply(pd.to_numeric).values + eps
            patient_data_feat_short_reads = patient_data[patient_data["Feature"]=="shortFragGCPara"][patient].apply(pd.to_numeric).values
            patient_data_feat_ultralong_reads = patient_data[patient_data["Feature"] == "ultraLongFragGCPara"][patient].apply(
                pd.to_numeric).values
            features_p.append([patient_data_feat_short_reads/patient_data_feat_reads,
                               patient_data_feat_ultralong_reads/patient_data_feat_reads
                               ])
        labels = metadata_patient[labels_cols].values
        wrapperDataTraining.append([np.asarray(features_p).T, labels])
        labels_training.append(labels)
        features_training.append(np.asarray(features_p).T)

    wrapperDataValidation = []
    labels_validation = []
    features_validation = []
    for patient in validation_columns:
        patient_data = data_feat[["Feature", patient]]
        metadata_patient = metadata[metadata["PatientId"] == patient]
        features_p = []
        for feat in features:
            # patient_data_feat = patient_data[patient_data["Feature"] == feat][patient].values
            patient_data_feat_reads = patient_data[patient_data["Feature"] == "noReadsGCParagonFrag"][patient].apply(pd.to_numeric).values + eps
            patient_data_feat_short_reads = patient_data[patient_data["Feature"] == "shortFragGCPara"][patient].apply(pd.to_numeric).values
            patient_data_feat_ultralong_reads = patient_data[patient_data["Feature"] == "ultraLongFragGCPara"][
                patient].apply(
                pd.to_numeric).values
            features_p.append([patient_data_feat_short_reads / patient_data_feat_reads,
                               patient_data_feat_ultralong_reads / patient_data_feat_reads
                               ])
        labels = metadata_patient[labels_cols].values
        wrapperDataValidation.append([np.asarray(features_p).T, labels])
        labels_validation.append(labels)
        features_validation.append(np.asarray(features_p).T)

    labels_training = np.asarray(labels_training).squeeze()
    features_training = np.asarray(features_training).squeeze()
    labels_validation = np.asarray(labels_validation).squeeze()
    features_validation = np.asarray(features_validation).squeeze()


    return wrapperDataTraining, wrapperDataValidation, labels_training, features_training, labels_validation,\
           features_validation, pt


file = "/Users/alexandra/PhD/FragmentationPatterns/Data/NN_5MB/NN_5MB_pivot_scaled.csv"
test_studies = ["SNR", "Normal"]
validation_studies = ["PRL", "Genome-IJB"]
features = ["noReadsGCParagonFrag", "shortFragGCPara"]

# get_data_wrapper(file, test_studies, validation_studies, features)

def get_data_wrapper_mixed(file, task, features):
    metadata = pd.read_csv("/Users/alexandra/PhD/FragmentationPatterns/Data/MetaData/AllMetaData.csv", sep="\t")
    labels_cols = ["PatientId", "Label", "ctDNADetected", "VAF"]
    data = pd.read_csv(file, sep="\t")
    data = data.fillna(0)
    metadata = metadata[(metadata["PatientId"].str.contains("PRL")) | (metadata["PatientId"].str.contains("SNR"))]

    # if task == "ctDNADetected" :
    #     check = 0
    # elif task =="Label":
    #     check="Healthy"
    # data_feat = data[data["Feature"].isin(features)]
    data_feat = data
    metadata = metadata[metadata[task].notna()]
    metadata["ctDNADetected"] = metadata["ctDNADetected"].apply(float)
    metadata["Label"] = metadata["Label"].apply(float)
    metadata["VAF"] = metadata["VAF"].apply(float)
    patients = list(metadata["PatientId"].values)
    unique_patients = u.get_list_of_patients(patients)
    metadata["UniquePatient"] = unique_patients
    healthy_patients = metadata[metadata[task] == 0]
    cancer_patients = metadata[metadata[task] == 1]

    healthy_metadata_train = healthy_patients.sample(round(len(healthy_patients)*0.6))
    cancer_metadat_train = cancer_patients.sample(round(len(cancer_patients)*0.6))
    train_dataset = pd.concat([healthy_metadata_train, cancer_metadat_train], axis=0)

    healthy_metadata_validation=healthy_patients[~healthy_patients["UniquePatient"].isin(train_dataset["UniquePatient"])]
    cancer_metadata_validation=cancer_patients[~cancer_patients["UniquePatient"].isin(train_dataset["UniquePatient"])]
    validation_dataset = pd.concat([healthy_metadata_validation, cancer_metadata_validation], axis=0)

    wrapperDataTraining, wrapperDataValidation, labels_training, features_training, labels_validation, \
    features_validation = create_data_wrapper(patients, data_feat, train_dataset["PatientId"].values, validation_dataset["PatientId"].values,
                metadata, labels_cols,features)

    return wrapperDataTraining, wrapperDataValidation, labels_training, features_training, labels_validation, features_validation

def get_data_wrapper_synergy_train(file, task, features, validation_study):
    metadata = pd.read_csv("/Users/alexandra/PhD/FragmentationPatterns/Data/MetaData/AllMetaData.csv", sep="\t")
    labels_cols = ["PatientId", "Label", "ctDNADetected", "VAF"]
    data = pd.read_csv(file, sep="\t")
    data = data.fillna(0)

    data_feat = data
    metadata = metadata[metadata[task].notna()]
    metadata["ctDNADetected"] = metadata["ctDNADetected"].apply(float)
    metadata["Label"] = metadata["Label"].apply(float)
    metadata["VAF"] = metadata["VAF"].apply(float)
    patients = list(metadata["PatientId"].values)

    train_dataset = metadata[((metadata["study"]=="Synergy") & (metadata["ctDNADetected"]==1)) | (metadata["PatientId"].str.contains("PIJB"))]

    validation_dataset = metadata[(metadata["study"].isin(validation_study)) | (metadata["PatientId"].str.contains("Genome-IJB"))]

    wrapperDataTraining, wrapperDataValidation, labels_training, features_training, labels_validation, features_validation\
        = create_data_wrapper(patients, data_feat, train_dataset["PatientId"].values,
                                              validation_dataset["PatientId"].values,
                                              metadata, labels_cols, features)

    return wrapperDataTraining, wrapperDataValidation, labels_training, features_training, labels_validation, features_validation


def get_data_wrapper_corrected(dagip_file, training_studies, validation_studies, task):
    metadata = pd.read_csv("/Users/alexandra/PhD/FragmentationPatterns/Data/MetaData/AllMetaData.csv", sep="\t")
    labels_cols = ["PatientId", "Label", "ctDNADetected", "VAF", "ichorTF"]
    data = pd.read_csv(dagip_file, sep="\t")
    data = data.fillna(0)

    # data_feat = data[data["Feature"].isin(features)]
    data_feat = data
    metadata = metadata[metadata[task].notna()]
    metadata["ctDNADetected"] = metadata["ctDNADetected"].apply(float)
    metadata["Label"] = metadata["Label"].apply(float)
    metadata["VAF"] = metadata["VAF"].apply(float)
    metadata["ichorTF"] = metadata["ichorTF"].apply(float)*100

    ichorTF = metadata["ichorTF"].values
    pt = PowerTransformer(method='yeo-johnson')
    ichorTF = ichorTF.reshape(-1, 1)
    pt.fit(ichorTF)
    transformed_data = pt.transform(ichorTF)
    metadata["ichorTFND"] = transformed_data

    # temporary fix see why this sample is missing
    patients = list(metadata[metadata["PatientId"] != "NIPT-PRL-062-25-D14-Px_S83"]["PatientId"].values)
    data_feat_cols = ["bin"] + ["Feature"] + patients
    data_feat = data_feat[data_feat_cols]
    all_columns = data_feat.columns
    eps = 0.001

    training_columns = []
    for training_study in training_studies:
        training_columns = training_columns + list(filter(lambda x: training_study in x, all_columns))

    validation_columns = []
    for validation_study in validation_studies:
        validation_columns = validation_columns + list(filter(lambda x: validation_study in x, all_columns))

    wrapperDataTraining = []
    labels_training = []
    features_training = []
    for patient in training_columns:
        features_p = []

        patient_data = data_feat[["Feature", patient]]
        metadata_patient = metadata[metadata["PatientId"] == patient]

        shortFragmentsRatio = patient_data[patient_data["Feature"] == "shortFragmentsRatio"][patient]
        longFragmentsRatio = patient_data[patient_data["Feature"] == "longFragmentsRatio"][patient]

        features_p.append([shortFragmentsRatio, longFragmentsRatio])

        labels = metadata_patient[labels_cols].values
        wrapperDataTraining.append([np.asarray(features_p).T, labels])
        labels_training.append(labels)
        features_training.append(np.asarray(features_p).T)

    wrapperDataValidation = []
    labels_validation = []
    features_validation = []
    for patient in validation_columns:
        features_p = []
        patient_data = data_feat[["Feature", patient]]
        metadata_patient = metadata[metadata["PatientId"] == patient]

        shortFragmentsRatio = patient_data[patient_data["Feature"] == "shortFragmentsRatio"][patient]
        longFragmentsRatio = patient_data[patient_data["Feature"] == "longFragmentsRatio"][patient]
        features_p.append([shortFragmentsRatio, longFragmentsRatio])

        labels = metadata_patient[labels_cols].values
        wrapperDataValidation.append([np.asarray(features_p).T, labels])
        labels_validation.append(labels)
        features_validation.append(np.asarray(features_p).T)

    labels_training = np.asarray(labels_training).squeeze()
    features_training = np.asarray(features_training).squeeze()
    labels_validation = np.asarray(labels_validation).squeeze()
    features_validation = np.asarray(features_validation).squeeze()

    return wrapperDataTraining, wrapperDataValidation, labels_training, features_training, labels_validation, \
           features_validation, pt


inp_file ="/Users/alexandra/PhD/PyCharmProjects/DAGIP/data/corrected_data_dagip.csv"



def get_train_syntethic_data(file, training_studies, validation_studies, features, task):
    # metadata = pd.read_csv("/Users/alexandra/PhD/FragmentationPatterns/Data/MetaData/AllMetaData.csv", sep="\t")
    labels_cols = ["PatientId", "Label", "ctDNADetected", "VAF"]
    no_reads = pd.read_csv("/Users/alexandra/PhD/FragmentationPatterns/Data/NN_5MB/NN_5MB_pivot_synth_noReadsGCParagonFrag_synergy.csv", sep="\t")
    short_reads = pd.read_csv(
        "/Users/alexandra/PhD/FragmentationPatterns/Data/NN_5MB/NN_5MB_pivot_synth_shortFragGCPara_synergy.csv",
        sep="\t")
    long_reads = pd.read_csv(
        "/Users/alexandra/PhD/FragmentationPatterns/Data/NN_5MB/NN_5MB_pivot_synth_ultraLongFragGCPara_synergy.csv",
        sep="\t")

    bins = list(range(1, 205, 1))


def get_nn_bias_corrected_data(task, features):

    test_data = pd.read_csv("/Users/alexandra/PhD/FragmentationPatterns/Data/NN_5MB/denosing_prediction/SRN_test_denoising_predicted.csv", sep="\t")
    validation_data = pd.read_csv("/Users/alexandra/PhD/FragmentationPatterns/Data/NN_5MB/denosing_prediction/PRL_validation_denoising_predicted.csv", sep="\t")

    metadata = pd.read_csv("/Users/alexandra/PhD/FragmentationPatterns/Data/MetaData/AllMetaData.csv", sep="\t")
    labels_cols = ["PatientId", "Label", "ctDNADetected", "VAF", "ichorTF"]
    data = pd.read_csv(file, sep="\t")
    data = data.fillna(0)

    # data_feat = data[data["Feature"].isin(features)]
    data_feat = data
    metadata = metadata[metadata[task].notna()]
    metadata["ctDNADetected"] = metadata["ctDNADetected"].apply(float)
    metadata["ichorTF"] = metadata["ichorTF"].apply(float)
    metadata["ichorTF"] = np.where(metadata["ichorTF"]<0.03, 0 ,metadata["ichorTF"]*100)
    metadata["Label"] = metadata["Label"].apply(float)
    metadata["VAF"] = metadata["VAF"].apply(float)

    ## transform nn output
    ichorTF = metadata["ctDNADetected"].values
    pt = PowerTransformer(method='yeo-johnson')
    ichorTF = ichorTF.reshape(-1, 1)
    pt.fit(ichorTF)
    transformed_data = pt.transform(ichorTF)
    metadata["ichorTFND"] = transformed_data

    # temporary fix see why this sample is missing
    missing = ["Genome-IJB-HP-40-xx_S28", "NIPT-PIJB-HP-73-Normal-15052023_S143"] #
    meta_v2 = metadata[~metadata["PatientId"].isin(missing)]
    patients = list(metadata[metadata["PatientId"] != "NIPT-PRL-062-25-D14-Px_S83"]["PatientId"].values)
    # patients = list(metadata["PatientId"].values)
    # patients = list(meta_v2["PatientId"].values)

    data_feat_cols = ["bin"] + ["Feature"] + patients
    data_feat = data_feat[data_feat_cols]
    all_columns = data_feat.columns
    eps = 0.001

    training_columns = test_data.columns
    validation_columns = validation_data.columns

    wrapperDataTraining = []
    labels_training = []
    features_training = []
    for patient in training_columns:
        patient_data = test_data[[patient]]
        metadata_patient = metadata[metadata["PatientId"] == patient]
        features_p = []
        for feat in features:
            features_p.append([patient_data[0:204],
                               patient_data[204:408]
                               ])
        labels = metadata_patient[labels_cols].values
        wrapperDataTraining.append([np.asarray(features_p).T, labels])
        labels_training.append(labels)
        features_training.append(np.asarray(features_p).T)

    wrapperDataValidation = []
    labels_validation = []
    features_validation = []
    for patient in validation_columns:
        patient_data = validation_data[[patient]]
        metadata_patient = metadata[metadata["PatientId"] == patient]
        features_p = []
        for feat in features:
            features_p.append([patient_data[0:204],
                               patient_data[204:408]
                               ])
        labels = metadata_patient[labels_cols].values
        wrapperDataValidation.append([np.asarray(features_p).T, labels])
        labels_validation.append(labels)
        features_validation.append(np.asarray(features_p).T)

    labels_training = np.asarray(labels_training).squeeze()
    features_training = np.asarray(features_training).squeeze()
    labels_validation = np.asarray(labels_validation).squeeze()
    features_validation = np.asarray(features_validation).squeeze()

    return wrapperDataTraining, wrapperDataValidation, labels_training, features_training, labels_validation, \
           features_validation, pt


def get_data_wrapper_2(file, task):
    metadata = pd.read_csv("/Users/alexandra/PhD/FragmentationPatterns/Data/MetaData/AllMetaData_.csv", sep="\t")
    labels_cols = ["PatientId", "Label", "ctDNADetected", "VAF", "ichorTF"]
    data = pd.read_csv(file, sep="\t")
    data = data.fillna(0)

    data_feat = data
    metadata = metadata[metadata[task].notna()]
    metadata["ctDNADetected"] = metadata["ctDNADetected"].apply(float)
    metadata["ichorTF"] = metadata["ichorTF"].apply(float)
    metadata["ichorTF"] = np.where(metadata["ichorTF"]<0.03, 0 ,metadata["ichorTF"]*100)
    metadata["Label"] = metadata["Label"].apply(float)
    metadata["VAF"] = metadata["VAF"].apply(float)

    ## transform nn output
    ichorTF = metadata["ichorTF"].values
    pt = PowerTransformer(method='yeo-johnson')
    ichorTF = ichorTF.reshape(-1, 1)
    pt.fit(ichorTF)
    transformed_data = pt.transform(ichorTF)
    metadata["ichorTFND"] = transformed_data
    patients = list(metadata["PatientId"].values)

    data_feat_cols = ["bin"] + ["Feature"] + patients
    data_feat = data_feat[data_feat_cols]

    training_columns = list(metadata[metadata["Process"]=="Training"]["PatientId"].values)

    validation_columns = list(metadata[metadata["Process"]=="Validation"]["PatientId"].values)

    wrapperDataTraining = []
    labels_training = []
    features_training = []
    for patient in training_columns:
        patient_data = data_feat[["Feature", patient]]
        metadata_patient = metadata[metadata["PatientId"]==patient]
        shortRatio = [patient_data[patient_data["Feature"] == "ratioShort"][patient].apply(pd.to_numeric).values]
        longRatio = [patient_data[patient_data["Feature"] == "ratioLong"][patient].apply(pd.to_numeric).values]
        labels = metadata_patient[labels_cols].values
        # wrapperDataTraining.append([np.asarray([shortRatio,longRatio]).T, labels])
        wrapperDataTraining.append([np.asarray(shortRatio).T, labels])
        labels_training.append(labels)
        # features_training.append(np.asarray([shortRatio,longRatio]).T)
        features_training.append(np.asarray(shortRatio).T)

    wrapperDataValidation = []
    labels_validation = []
    features_validation = []
    for patient in validation_columns:
        patient_data = data_feat[["Feature", patient]]
        metadata_patient = metadata[metadata["PatientId"] == patient]
        shortRatio = [patient_data[patient_data["Feature"] == "ratioShort"][patient].apply(pd.to_numeric).values]
        longRatio = [patient_data[patient_data["Feature"] == "ratioLong"][patient].apply(pd.to_numeric).values]
        labels = metadata_patient[labels_cols].values
        # wrapperDataValidation.append([np.asarray([shortRatio,longRatio]).T, labels])
        wrapperDataValidation.append([np.asarray(shortRatio).T, labels])
        labels_validation.append(labels)
        # features_validation.append(np.asarray([shortRatio,longRatio]).T)
        features_validation.append(np.asarray(shortRatio).T)

    labels_training = np.asarray(labels_training).squeeze()
    features_training = np.asarray(features_training).squeeze()
    labels_validation = np.asarray(labels_validation).squeeze()
    features_validation = np.asarray(features_validation).squeeze()


    return wrapperDataTraining, wrapperDataValidation, labels_training, features_training, labels_validation,\
           features_validation, pt

def get_data_insilico(file, task):
    metadata = pd.read_csv("/Users/alexandra/PhD/PyCharmProjects/ALFAssay/validations/data/synthetic_AllMetaData.csv", sep="\t")
    labels_cols = ["PatientId", "Label", "ctDNADetected", "VAF", "ichorTF"]
    data = pd.read_csv(file, sep="\t")
    data = data.fillna(0)

    data_feat = data
    metadata = metadata[metadata[task].notna()]
    metadata["ctDNADetected"] = metadata["ctDNADetected"].apply(float)
    metadata["ichorTF"] = metadata["ichorTF"].apply(float)
    metadata["ichorTF"] = np.where(metadata["ichorTF"] < 0.03, 0, metadata["ichorTF"] * 100)
    metadata["Label"] = metadata["Label"].apply(float)
    metadata["VAF"] = metadata["VAF"].apply(float)

    patients = list(metadata["PatientId"].values)

    data_feat_cols =  patients
    data_feat = data_feat[data_feat_cols]


    validation_columns = list(metadata[metadata["Process"] == "Validation"]["PatientId"].values)

    wrapperDataValidation = []
    labels_validation = []
    features_validation = []
    for patient in validation_columns:
        patient_data = data_feat[[patient]]
        metadata_patient = metadata[metadata["PatientId"] == patient]
        shortRatio = [patient_data.apply(pd.to_numeric).values]
        labels = metadata_patient[labels_cols].values
        # wrapperDataValidation.append([np.asarray([shortRatio,longRatio]).T, labels])
        wrapperDataValidation.append([np.asarray(shortRatio).T, labels])
        labels_validation.append(labels)
        # features_validation.append(np.asarray([shortRatio,longRatio]).T)
        features_validation.append(np.asarray(shortRatio).T)

    labels_validation = np.asarray(labels_validation).squeeze()
    features_validation = np.asarray(features_validation).squeeze()

    return  wrapperDataValidation, labels_validation, \
           features_validation