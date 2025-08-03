# Created by alexandra at 20/12/2023
from math import exp
import torch as t
import os.path
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler, MinMaxScaler

def sigmoid(x):
    return 1.0 / (1.0 + exp(-x))


def accuracy_fn(y_true, y_pred):
    correct = t.eq(y_true, y_pred).sum().item() # torch.eq() calculates where two tensors are equal
    acc = (correct / len(y_pred)) * 100
    return acc

def get_patient_from_file_name(filename, pattern):
    patient_id = os.path.basename(filename).replace(pattern, "")
    return patient_id

def basic_scale(x):
    return (x - np.mean(x))/np.std(x)


def norm_alternative(x):
    return (1e6*x/np.sum(x))

def norm_min_max(x):
    return (x - np.min(x))/(np.max(x) - np.min(x))


def norm_up_one(x):
    return 1e6 * x/np.sum(x)


def getSampleIdSynthetic(sampleId):
    data = pd.read_csv("/Users/alexandra/PhD/FragmentationPatterns/Data/NN_5MB/NN_5MB_scaled.csv", sep="\t")
    sample_data = data[data["PatientId"]==sampleId]
    sampleBatch = sample_data["batch"].values
    fragmentSize = []
    fragmentSize += [t.from_numpy(np.array(sample_data["noReadsGCParagonFrag", "shortOverLongGCPara"].values, dtype=np.float64))]
    fragmentSize = t.stack(fragmentSize)
    sampleBatches = t.stack(sampleBatch)

    return fragmentSize, sampleBatches, sampleId

#
def scale(dataset):
    scaler = MinMaxScaler()

    return scaler


def fix_clustermap(plot, title):
    plot.fig.suptitle(title)
    # plot.cax.set_visible(False)
    plt.subplots_adjust(bottom=0.1)
    # Remove the dendrogram (https://stackoverflow.com/a/42418968/1878788)
    plot.ax_row_dendrogram.set_visible(False)
    # Use the dendrogram box to reposition the colour bar
    dendro_box = plot.ax_row_dendrogram.get_position()
    dendro_box.x0 = (dendro_box.x0 + 2 * dendro_box.x1) / 3
    plot.cax.set_position(dendro_box)
    # Move the ticks to the left (https://stackoverflow.com/a/36939552/1878788)
    plot.cax.yaxis.set_ticks_position("left")

    hm = plot.ax_heatmap.get_position()
    plt.setp(plot.ax_heatmap.yaxis.get_majorticklabels(), fontsize=5)
    plot.ax_heatmap.set_position([hm.x0, hm.y0, hm.width, hm.height])


def get_patient_from_sample(sample):
    if "PRL" in sample: #NIPT-PRL-001-25-BL-Px_S1
        patient = str(sample.strip('NIPT-PRL-').split("-")[0].split("_")[0]) + "PRL"
    elif "NRH" in sample: #NIPT-NRH-027-PL3-Nx_S72
        patient = str(sample.strip('NIPT-NRH-').split("-")[0].split("_")[0]) + "NRH"
    elif "SNR" in sample: #NIPT-P0001-ArmPhaseI-SNRBE0025-W1-Screening_S118
        patient = str(sample.strip('NIPT-').split("-")[0].split("_")[0]) + "SNR"
    elif "PIJB" in sample:
        patient = str(sample.strip('NIPT-PIJB-HP-').split("-")[0].split("_")[0]) + "SH"
    elif "Genome" in sample:
        patient = str(sample.strip('Genome-IJB-HP-').split("-")[0].split("_")[0]) + "H"
    else:
        patient=""
    return patient

def get_list_of_patients(samples):
    patients = []
    for sample in samples:
        patient = get_patient_from_sample(sample)
        patients.append(patient)
    return patients


def RMSELoss(yhat, y):
    return t.sqrt(t.mean((yhat-y)**2))