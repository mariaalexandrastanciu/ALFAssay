# Created by alexandra at 01/12/2022
import numpy as np
import glob
import pandas as pd
from FeatureEngineering import CONSTANTS as C, utils as u
import os.path
from pathlib import Path


class DataWrapper:
    def __init__(self, cancer_paths, healthy_path, output_path, gather_samples=False, ichorCNAThreshold=0,
                 VAFThreshold=0.2, start_fragment_count_size=30,
                 start_fragment=30, end_fragment=700, genome_version="hg19", task=C.Label):
        self.cancer_paths = cancer_paths
        self.healthy_path = healthy_path + "/*.csv"
        self.output_path = output_path
        self.gather_samples = gather_samples
        self.ichorCNAThreshold = ichorCNAThreshold
        self.VAFThreshold = VAFThreshold
        studies = []
        output_gather_frag_size_by_study = []
        output_gather_frag_size_by_ichor = []
        output_gather_frag_size_by_vaf = []
        for cancer_path in cancer_paths:
            study = u.get_study_from_path(cancer_path)
            output_gather_frag_size_by_study.append(self.output_path + "/" + study + "_" + C.fragSizeNaming)
            output_gather_frag_size_by_ichor.append(self.output_path + "/" + study + "_ichor_" + C.fragSizeNaming)
            output_gather_frag_size_by_vaf.append(self.output_path + "/" + study + "_vaf_" + C.fragSizeNaming)
            studies.append(study)
        self.studies = studies
        self.output_gather_frag_size_by_studies = output_gather_frag_size_by_study
        self.output_gather_frag_size_by_ichor = output_gather_frag_size_by_ichor
        self.output_gather_frag_size_by_vaf = output_gather_frag_size_by_vaf
        self.start_fragment = start_fragment
        self.end_fragment = end_fragment
        self.start_fragment_count_size = start_fragment_count_size
        self.genome_version = genome_version
        if task == C.ctDNADetected:
            self.field = "ctDNADetected"
            self.value_true = "Detected"
            self.value_false = "NotDetected"
        elif task == C.VAFg0p001:
            self.field = "VAFg0p001"
            self.value_true = "Detected"
            self.value_false = "NotDetected"
        else:
            self.field = "Label"
            self.value_true = "Cancer"
            self.value_false = "Healthy"

    def gather_data_for_analysis(self):
        if self.gather_samples:
            self.gather_frag_size_samples_by_study()
            self.group_by_all_timepoints()
            self.group_by_ichor_threshold()
            self.group_by_VAF_threshold()
            self.group_by_ctDNA_detection()
            self.group_by_TumorAssesment()
            # self.create_delfi_data()
            self.group_by_VAFTumorAssesment()

    def group_by_all_timepoints(self):
        for i, cancer_path in enumerate(self.cancer_paths):
            frag_size = pd.read_csv(self.output_gather_frag_size_by_studies[i] + ".csv", sep="\t")
            for timepoint in C.timepoints:
                frag_size_by_timepoint = frag_size[(frag_size["timepoint"] == timepoint) |
                                                   (frag_size["Label"] == "Healthy")]
                frag_size_by_timepoint.to_csv(self.output_gather_frag_size_by_studies[i] + "_" + timepoint + ".csv",
                                              sep="\t", index=False, header=True)

    def group_by_ichor_threshold(self):
        for i, cancer_path in enumerate(self.cancer_paths):
            frag_size = pd.read_csv(self.output_gather_frag_size_by_studies[i] + ".csv", sep="\t")
            frag_size_by_ichor = frag_size[(frag_size["ichorTF"] >= self.ichorCNAThreshold) |
                                           (frag_size["Label"] == "Healthy")]
            frag_size_by_ichor.to_csv(self.output_gather_frag_size_by_ichor[i] + "_" + str(self.ichorCNAThreshold) +
                                      ".csv", sep="\t", index=False, header=True)

    def group_by_VAF_threshold(self):
        for i, cancer_path in enumerate(self.cancer_paths):
            frag_size = pd.read_csv(self.output_gather_frag_size_by_studies[i] + ".csv", sep="\t")
            frag_size_by_VAF = frag_size[(frag_size["VAF"] >= self.VAFThreshold) |
                                         (frag_size["Label"] == "Healthy")]
            frag_size_by_VAF.to_csv(self.output_gather_frag_size_by_vaf[i] + "_" + str(self.VAFThreshold) + ".csv",
                                    sep="\t", index=False, header=True)

    def group_by_ctDNA_detection(self):
        for i, cancer_path in enumerate(self.cancer_paths):
            frag_size = pd.read_csv(self.output_gather_frag_size_by_studies[i] + ".csv", sep="\t")
            frag_size_by_ctDNADetection = frag_size[(frag_size["ctDNADetected"] == True) |
                                                    (frag_size["Label"] == "Healthy")]
            frag_size_by_ctDNADetection.to_csv(self.output_gather_frag_size_by_vaf[i] + "_" + "ctDNATrue" + ".csv",
                                               sep="\t", index=False, header=True)


    def group_by_TumorAssesment(self):
        for i, cancer_path in enumerate(self.cancer_paths):
            frag_size = pd.read_csv(self.output_gather_frag_size_by_studies[i] + ".csv", sep="\t")
            frag_size_by_ctDNADetection = frag_size[frag_size['ctDNADetected'].notna()]
            frag_size_by_ctDNADetection['ctDNADetected'] = np.where(frag_size_by_ctDNADetection['ctDNADetected'] == True,
                                                                    "Detected", "NotDetected")

            frag_size_by_ctDNADetection.to_csv(self.output_gather_frag_size_by_vaf[i] + "_" + "ctDNADetection" + ".csv",
                                               sep="\t", index=False, header=True)

    def group_by_VAFTumorAssesment(self):
        for i, cancer_path in enumerate(self.cancer_paths):
            frag_size = pd.read_csv(self.output_gather_frag_size_by_studies[i] + ".csv", sep="\t")

            frag_size_by_ctDNADetection = frag_size[frag_size['VAFg0p001'].notna()]
            frag_size_by_ctDNADetection.to_csv(self.output_gather_frag_size_by_vaf[i] + "_" + "VAFg0p001" + ".csv",
                                               sep="\t", index=False, header=True)

    def gather_frag_size_samples_by_study(self):
        for i, cancer_path in enumerate(self.cancer_paths):
            cols = []
            frag_size_cancer_files_all = glob.glob(cancer_path + C.fragmentSizeFolder + "/*.csv")

            new_cols_frags = list(range(self.start_fragment, self.end_fragment+1))
            cols_frags = list(range(self.start_fragment - self.start_fragment_count_size,
                                    self.end_fragment - self.start_fragment_count_size+1))

            h_samples_df = self.add_healthy(new_cols_frags, cols_frags)
            cancer_samples = []

            for frag_size_cancer_file_all in frag_size_cancer_files_all:
                patient_id = u.get_patient_from_file_name(frag_size_cancer_file_all,
                                                          "_fragment_size_summary.csv")


                all_frag_size = pd.read_csv(frag_size_cancer_file_all, sep="\t")[["fragment_counts"]].T
                frag_size = all_frag_size[cols_frags]
                frag_size.columns = new_cols_frags

                meta_data = pd.read_csv("/Users/alexandra/PhD/FragmentationPatterns/Data/MetaData/AllMetaDataLabels.csv", sep="\t")
                meta_data_patient = meta_data[meta_data["PatientId"]==patient_id]

                frag_size["PatientId"] = patient_id
                frag_size["Label"] = meta_data_patient["Label"].values
                frag_size["ichorTF"] = meta_data_patient["ichorTF"].values
                frag_size["timepoint"] = meta_data_patient["timepoint"].values
                frag_size["VAF"] = meta_data_patient["VAF"].values
                frag_size["ctDNADetected"] = meta_data_patient["ctDNADetected"].values
                frag_size["survivalStatus"] = meta_data_patient["survivalStatus"].values
                frag_size["survivalTime"] = meta_data_patient["survivalTime"].values
                frag_size["VAFg0p001"] = meta_data_patient["VAFg0p001"].values
                frag_size["study"] = meta_data_patient["study"].values
                # frag_size["batch"] = meta_data_patient["batch"].values
                cancer_samples.append(frag_size.values)
                cols = frag_size.columns

            cancer_samples_flat = [item for sublist in cancer_samples for item in sublist]
            cancer_samples_df = pd.DataFrame(cancer_samples_flat, columns=cols)
            all_samples_df = pd.concat([cancer_samples_df, h_samples_df])
            all_samples_df.to_csv(self.output_gather_frag_size_by_studies[i] + ".csv", sep="\t",
                                  index=False, header=True)

    def add_healthy(self, new_cols_frags, cols_frags):
        frag_size_healthy_files = glob.glob(self.healthy_path)
        cols = []
        h_samples = []
        for frag_size_healthy_file in frag_size_healthy_files:
            patient_id = u.get_patient_from_file_name(frag_size_healthy_file, "_fragment_size_summary.csv")
            frag_size = pd.read_csv(frag_size_healthy_file, sep="\t")[["fragment_counts"]].T
            frag_size = frag_size[cols_frags]
            frag_size.columns = new_cols_frags

            frag_size["PatientId"] = patient_id
            frag_size["Label"] = "Healthy"
            frag_size["ichorTF"] = 0
            frag_size["timepoint"] = "NA"
            frag_size["VAF"] = 0
            frag_size["ctDNADetected"] = "NotDetected"
            frag_size["survivalStatus"] = 0
            frag_size["survivalTime"] = "NA"
            frag_size["VAFg0p001"] = "NotDetected"
            frag_size["study"] = "Healthy"
            # frag_size["batch"] = meta_data_patient["batch"].values

            cols = frag_size.columns
            h_samples.append(frag_size.values)
        h_samples_flat = [item for sublist in h_samples for item in sublist]
        h_samples_df = pd.DataFrame(h_samples_flat, columns=cols)

        return h_samples_df


    def filter_fragment_lengths(self, full_fragment_length_path, initial_start_fragment_length,
                                initial_end_fragment_length, filter_start_fragment_length,
                                filter_end_fragment_length, outputdir):

        full_fragment_length_files = glob.glob(full_fragment_length_path + '/*csv')
        start = filter_start_fragment_length - initial_start_fragment_length
        end = start + filter_end_fragment_length - filter_start_fragment_length
        for full_fragment_length_file in full_fragment_length_files:
            file_name = os.path.basename(full_fragment_length_file)
            data = pd.read_csv(full_fragment_length_file, sep="\t")
            data_list = data.values

            features_filtered = data_list[:, start:end]
            labels = data_list[:, -8:]
            cols = data.columns[0:end - start].tolist() + data.columns[-8:].tolist()
            data_filtered_array = np.concatenate([features_filtered, labels], axis=1)
            data_filtered = pd.DataFrame(data_filtered_array, columns=cols)
            output_filename = outputdir + "/" + file_name
            data_filtered.to_csv(output_filename, sep="\t", index=False, header=True)


    def create_delfi_data(self):
        delfi_data = []
        meta_data_delfi = pd.read_csv("/Users/alexandra/PhD/FragmentationPatterns/Data/MetaData/DELFIMetaData.csv",
                                      sep="\t")
        patient_list = meta_data_delfi["PatientId"].values
        new_cols_frags = list(range(self.start_fragment, self.end_fragment + 1))
        cols_frags = list(range(self.start_fragment - self.start_fragment_count_size,
                                self.end_fragment - self.start_fragment_count_size + 1))

        for patientId in patient_list:
            meta_data_delfi_per_patients = meta_data_delfi[meta_data_delfi["PatientId"] == patientId]
            file = "/Users/alexandra/PhD/DELFIStudy/FragmentationPatterns/FragSize_30_700_mq40/" + \
                   patientId + "_fragment_size_summary.csv"
            my_file = Path(file)
            if my_file.is_file():
                fragmentSizes_pd = pd.read_csv(file, sep="\t")[["fragment_counts"]].T
                fragmentSizes = fragmentSizes_pd[cols_frags]
                fragmentSizes.columns = new_cols_frags
                # fragmentSizes[meta_data_delfi_per_patients.columns] = meta_data_delfi_per_patients
                fragmentSizes["PatientId"] = patientId
                fragmentSizes["Label"] = meta_data_delfi_per_patients["Label"].values
                fragmentSizes["ichorTF"] = "NA"
                fragmentSizes["timepoint"] = "BL"
                fragmentSizes["VAF"] = "NA"
                fragmentSizes["ctDNADetected"] = "NA"
                fragmentSizes["survivalStatus"] = "NA"
                fragmentSizes["survivalTime"] = "NA"
                delfi_data.append(fragmentSizes.values)


        delfi_data_flat = [item for sublist in delfi_data for item in sublist]
        all_cols = new_cols_frags + meta_data_delfi.columns.tolist()
        delfi_data_df = pd.DataFrame(delfi_data_flat, columns=all_cols)
        delfi_data_df.to_csv(self.output_path + "/" + "DELFIStudy" + "_" + C.fragSizeNaming + ".csv", sep="\t",
                             index=False, header=True)
        print("ok")

