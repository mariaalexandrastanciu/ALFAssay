# Created by alexandra at 19/03/2025
import pandas as pd
import os.path
import Visualisation as V
import numpy as np


def create_predictions_with_extra_info(pred_file, meta):
    predictions = pd.read_csv(pred_file, sep="\t")
    metadata = pd.read_csv(meta, sep="\t")

    pred_with_extra_info = pd.merge(predictions, metadata, on="PatientId", how="inner")
    file = os.path.basename(pred_file).strip('.csv')
    pred_with_extra_info.to_csv("/Users/alexandra/PhD/PyCharmProjects/ALFAssayNN/results/" + file + ".csv" , sep="\t", index=False)

# pred_file="/Users/alexandra/PhD/PyCharmProjects/ALFAssayNN/plots/results_validation_PRL.csv"
# meta="/Users/alexandra/PhD/FragmentationPatterns/Data/MetaData/AllMetaData.csv"
# create_predictions_with_extra_info(pred_file, meta)


def roc_with_ichor(pred_file, ichorCNA_file):
    predictions = pd.read_csv(pred_file, sep="\t")
    ichorCNA = pd.read_csv(ichorCNA_file, sep="\t")
    # predictions = predictions[(predictions["PatientId"].str.contains("NRH")) | (predictions["PatientId"].str.contains("IJB"))]
    predictions = predictions[(predictions["PatientId"].str.contains("SNR"))]


    data = pd.merge(predictions, ichorCNA, on="PatientId")


    # ichorCNA = ichorCNA[ichorCNA["SampleID"].isin(predictions["PatientId"])]

    title="ALFAssay on ctDNA and ichorCNA performance on SNR data"
    output_file = "/Users/alexandra/PhD/PyCharmProjects/ALFAssayNN/plots/plot_prl_train/"
    caption = '''ALFAssay on ctDNADetection and ichorCNA performance on SNR data,\n model trained on PRL and healthy neorhea '''
    # caption = '''ALFAssay on ctDNAdetection and ichorCNA performance on PRL data,\n model trained on all SYN and healthy Synergy'''
    V.plot_roc_fpr_tpr_ichor(data["Truelabel"], data["Predictedlabel"], output_file, title, 1, data["ichorTF"], caption)




def check_pred_ichor(pred_file, ichorCNA_file):
    predictions = pd.read_csv(pred_file, sep="\t")
    ichorCNA = pd.read_csv(ichorCNA_file, sep="\t")
    # predictions = predictions[(predictions["PatientId"].str.contains("NRH")) | (predictions["PatientId"].str.contains("IJB"))]
    predictions = predictions[(predictions["PatientId"].str.contains("PRL"))]


    data = pd.merge(predictions, ichorCNA, on="PatientId")

    title = "ALFAssay on Survival status and ichorCNA performance on PRL data"
    output_file = "../NNModel/plots/SurvivalModelPlots/"
    caption = '''ALFAssay on Survival status and ichorCNA performance on PRL data,\n model trained on SNR cancer '''
    # caption = '''ALFAssay on ctDNAdetection and ichorCNA performance on PRL data,\n model trained on all SYN and healthy Synergy'''

    # true_y = np.where(data["ichorTF"] > 0.03, 1, 0)
    # pred_y = np.where(data["Predictedlabel"] > 0.03, 1, 0)
    V.plot_roc_fpr_tpr_ichor(data["Truelabel"], data["Predictedlabel"], output_file, title, 1, data["ichorTF"], caption)

pred_file="../NNModel/results/results_validation_NRH_PRL_healthy_survLabel.csv"
metadata="/Users/alexandra/PhD/FragmentationPatterns/Data/MetaData/AllMetaData.csv"
check_pred_ichor(pred_file, metadata)