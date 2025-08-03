# Created by alexandra at 21/11/2022
(PatientId, Label, ichorTF, timepoint, VAF, ctDNADetected, status, OSF, VAFg0p001) = range(9)
# chromosomes = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
#                "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22"]

chromosomes_dict = {"chr1":0, "chr2":1, "chr3":2, "chr4":3, "chr5":4, "chr6":5, "chr7":6, "chr8":7, "chr9":8, "chr10":9,
                    "chr11":10, "chr12":11, "chr13":12, "chr14":13, "chr15":14, "chr16":15, "chr17":16, "chr18":17,
                    "chr19":18, "chr20":19, "chr21":20, "chr22":21}
chromosomes = ["chr19", "chr21"]
#  removed in DELFI: 13p,14p,15p,21p,22p
# chr_arms = ["1p", "1q", "2p", "2q", "3p", "3q", "4p", "4q", "5p", "5q", "6p", "6q",
#             "7p", "7q", "8p", "8q", "9p", "9q", "10p", "10q", "11p", "11q", "12p",
#             "12q", "13p", "13q", "14p", "14q", "15p", "15q", "16p", "16q", "17p", "17q", "18p", "18q",
#              "19p", "19q", "20p", "20q", "21p", "21q", "22p", "22q"]
chr_arms = ["19p", "19q", "21p", "21q"]

non_somatic_chr = ["chrX", "chrM", "chrY"]

timepoints = ["BL", "C1", "Surgery"]

data_classes = {"healthy": 0, "cancer": 1}

fragmentPatternsFolder = "/FragmentationPatterns"
fragmentSizeMutationsFolder = "/FragmentationPatterns/FragmentSizes/FragSizeMutations"
# fragmentSizeFolder = "/FragmentationPatterns/FragmentSizes/FragSizeGCCorrected_hg38"
fragmentSizeFolder = "/FragmentationPatterns/FragmentSizes/NN_5MB_summary"
# fragmentSizeFolderHealthy = "/FragmentationPatterns/FragmentSizes/FragSizeGCCorrected_hg38"
fragmentSizeFolderHealthy = "/FragmentationPatterns/FragmentSizes/NN_5MB_summary"
fragmentSizeWindowFolder = "/FragmentationPatterns/WindowFragmentSizes"
#label_types
(sampleID, cancerType, ichorCNA) = range(3)

fragSizeNaming = "fragCounts"
fragSizeByIchorNaming = "fragCounts_by_ichor_"
fragSizeByTimepointNaming = "fragCounts_by_timepoint_"
fragSizeByStudyNaming = "fragCounts_by_study_"
analysis_types = ["byIchorCNA", "byStudy", "byTimepoint"]
fragCountsGatheredNaming = {"byIchorCNA": fragSizeByIchorNaming, "byTimepoint": fragSizeByTimepointNaming,
                            "byStudy": fragSizeByStudyNaming}
fragSizeParamUsage = {"byIchorCNA": "ichorThreshold", "byTimepoint": "timepoint", "byStudy": "study" }


task_naming = {Label: "Label", ctDNADetected: "ctDNADetected", VAFg0p001: "VAFg0p001"}
validationDatasetNeorhea = "Neorhea"
trainNeorheaValidatePearl = "Pearl"
validationDatasetMix = "mix"
validationDatasetMixTestValidNeorhea = "validNeorhea"
testAndValidateNeorhea = "Train_ValidNeorhea"
testPearlValidDELFI = "testPearlValidDELFI"
testNeorheaValidDELFI = "testNeorheaValidDELFI"
testNeorheaHealthyDelfiValidDELFIHealthyNeo = "testNeorheaHealtyDelfiValidDELFIHealthyNeo"
testNeorheaMixHealthyValidDELFIMixHealthy = "testNeorheaMixHealthyValidDELFIMixHealthy"
trainPearlValidateNeorheaHealthy="trainPearlValidateNeorheaHealthy"
trainNeorheaValidatePearlHealthy="trainHealthyValidatePearlHealthy"
testAndValidateDELFI = "ValidDELFI"
input_type_per_bin = "per_bin"
input_type_per_sample = "per_sample"
input = {"per_bin": "/window_fragsize_short_noReads_by_bin.csv", "per_sample": "/window_fragsize_short_noReads.csv"}

output_type = {"table": "Tables", "plot" : "Plots"}

task_values = {ctDNADetected: ["NotDetected", "Detected"], Label: ["Healthy", "Cancer"],
               VAFg0p001:["NotDetected", "Detected"]}
