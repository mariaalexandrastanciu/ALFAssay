# Created by alexandra at 28/11/2022
import numpy as np
import glob
import os.path
import os
import pandas as pd


def NormalizeData(data):
    return (data - np.min(data)) / (np.max(data) - np.min(data))


def getTimepoint_structure(file_name):
    timepoints = []

    if "NRH" in file_name or "NeoRhea" in file_name:
        timepoints = ["PL1", "PL2", "PL3"]
    elif 'PRL' in file_name or "Pearl" in file_name:
        timepoints = ["BL", "D14", "PD"]
    else:
        timepoints = [""]
    return timepoints


def get_patient_from_file_name(filename, pattern):
    patient_id = os.path.basename(filename).replace(pattern, "")
    return patient_id


def basic_scale(x):
    return (x - np.mean(x))/np.std(x)


def get_study_from_path(path):
    study = os.path.basename(path)
    return study

def get_timepoint_from_path(string):
    timepoint = "None"
    if "BL" in string or "PL1" in string or "W1" in string or "w1" in string:
        timepoint = "BL"
    if "D14" in string or "PL2" in string or "C1" in string or "W3" in string or "w3" in string:
        timepoint = "C1"
    if "PD" in string or "PL3" in string or "Surgery" in string or "W13" in string or "w13" in string:
        timepoint = "Surgery"
    return timepoint

def get_title_from_filename(filename):
    timepoint = get_timepoint_from_path(filename)
    if timepoint == "None":
        timepoint = ""
    if "NRH" in filename or "NeoRhea" in filename:
        study = "Neorhea"
    elif 'PRL' in filename or "Pearl" in filename:
        study = "Pearl"
    elif 'SRN' in filename or "Synergy" in filename:
        study = "Synergy"
    else:
        study = ""
    title = timepoint + " " + study
    if "ichor" in filename:
        title = timepoint + " " + study + " for ichorCNA with TF > 0.1"

    if "ctDNATrue" in filename and "vaf" in filename:
        title = timepoint + " " + study + " for ctDNA detected samples"
    elif "vaf" in filename:
        title = timepoint + " " + study + " for VAF values > 0.2"

    return title

def gather_frag_size_per_sample(input_files, output_file, start_fragment_size, end_fragment_size):
    gathered_data = []
    for file in input_files:
        sample_data = pd.read_csv(file, sep="\t")["fragment_counts"].values.T
        gathered_data.append(sample_data)


    columns = [str(x) for x in list(range(start_fragment_size, end_fragment_size+1))]
    gathered_data_df = pd.DataFrame(gathered_data, columns=columns)
    gathered_data_df.to_csv(output_file, sep="\t", index=False)


def get_label_healthy(string):
    label = "Healthy"
    if "PIJB" in string:
        label = "Healthy_Synergy"
    return label

def get_ctDNA_synergy(string):
    positive=0
    week1_ctDNA_pos = pd.read_csv("/Users/alexandra/PhD/SynergyStudy/metadata/week1_ctDNA_David.csv", sep=",")
    positive_samples = week1_ctDNA_pos[week1_ctDNA_pos["ctDNA."]=="Y"]["Sample"].values
    for positive_sample in positive_samples:
        if positive_sample in string:
            positive=1
    return positive