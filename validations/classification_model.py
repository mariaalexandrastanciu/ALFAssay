# Created by alexandra at 31/03/2025

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm

from sklearn import datasets
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process.kernels import RBF
from sklearn.inspection import DecisionBoundaryDisplay
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score
from sklearn.multiclass import OneVsRestClassifier
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
import pandas as pd
from sklearn.model_selection import RepeatedStratifiedKFold
import Visualisation as V


def perform_classification(test_snr, validation_prl, validation_nrh, x_feature_test, y_feature_test,
                           x_feature_validation, y_feature_validation, model):

    print("Accuracy for model: ", model)
    if len(x_feature_test) > 1 :
        X = test_snr[x_feature_test].values
        X_validation_prl = validation_prl[x_feature_validation].values
        X_validation_nrh = validation_nrh[x_feature_validation].values
        y = np.asarray(test_snr[y_feature_test].values, dtype=int).squeeze()  # test_snr[y_feature_test].apply(int).values
        y_validation_prl = np.asarray(validation_prl[y_feature_validation].values,
                                      dtype=int).squeeze() # data_validation_prl_all[y_feature_validation].apply(int).values
        y_validation_nrh = np.asarray(validation_nrh[y_feature_validation].values,
                                      dtype=int).squeeze()  # data_validation_nrh_all[y_feature_validation].apply(int).values
    else:
        X = test_snr[x_feature_test].values.reshape(-1, 1)
        X_validation_prl = validation_prl[x_feature_validation].values.reshape(-1, 1)
        X_validation_nrh = validation_nrh[x_feature_validation].values.reshape(-1, 1)
        y = np.asarray(test_snr[y_feature_test].values, dtype=int).squeeze() # test_snr[y_feature_test].apply(int).values
        y_validation_prl = np.asarray(validation_prl[y_feature_validation].values,
                                      dtype=int).squeeze()   # data_validation_prl_all[y_feature_validation].apply(int).values
        y_validation_nrh = np.asarray(validation_nrh[y_feature_validation].values, dtype=int).squeeze()



    # n_features = 1

    rskf = RepeatedStratifiedKFold(n_splits=5, n_repeats=1, random_state=36851234)
    C = 10
    # clf = SVC(kernel="linear", C=C, probability=True, random_state=0)
    clf = RandomForestClassifier(random_state=0)

    probabilities_test = []
    y_srn = []
    accuracy_cv = []
    for i, (train_index, test_index) in enumerate(rskf.split(X, y)):
        clf.fit(X[train_index], y[train_index])
        prob_test = clf.predict_proba(X[test_index])
        probabilities_test.extend(prob_test)
        y_pred_test = clf.predict(X[test_index])
        accuracy_test = accuracy_score(y[test_index], y_pred_test)
        y_srn.extend(y[test_index])
        accuracy_cv.extend([accuracy_test])

    print(f"Accuracy (test) for Synergy: {np.mean(accuracy_cv):0.1%}")

    probabilities_prl = clf.predict_proba(X_validation_prl)
    y_pred_prl = clf.predict(X_validation_prl)
    accuracy_prl = accuracy_score(y_validation_prl, y_pred_prl)
    print(f"Accuracy (validation) for Pearl: {accuracy_prl:0.1%}")

    probabilities_nrh = clf.predict_proba(X_validation_nrh)
    y_pred_nrh = clf.predict(X_validation_nrh)
    accuracy_nrh = accuracy_score(y_validation_nrh, y_pred_nrh)
    print(f"Accuracy (validation) for Neorhea: {accuracy_nrh:0.1%}")

    return probabilities_test, probabilities_prl, probabilities_nrh, y_srn, y_validation_prl, y_validation_nrh

pred_file_test = "../NNModel/results/results_train_SNR_healthy_ichorTF_t3.csv"
pred_file_validation ="../NNModel/results/results_validation_NRH_PRL_healthy_ichorTF_t3.csv"
metadata = pd.read_csv("/Users/alexandra/PhD/FragmentationPatterns/Data/MetaData/AllMetaData.csv", sep="\t")
ichor_healthy = pd.read_csv("/Users/alexandra/PhD/healthy_sWGS/ichor_CNA_results.csv", sep="\t")
metadata_healthy = metadata[metadata["Label"]==0].drop("ichorTF", axis=1)
metadata_cancer = metadata[metadata["Label"]!=0]
metadata_healthy = pd.merge(metadata_healthy, ichor_healthy, on="PatientId")[metadata_cancer.columns]
metadata = pd.concat([metadata_cancer, metadata_healthy], axis=0)


### fragle results ###
pred_file_fragle_prl = pd.read_csv("/Users/alexandra/PhD/PearlStudy/FragmentationPatterns/fragle/Fragle_Pearl.csv", sep=",")
# pred_file_fragle_prl = pd.concat([pred_file_fragle_prl, fragle_valid], axis=0)
pred_file_fragle_prl["PatientId"] = pred_file_fragle_prl["Sample_ID"].str.replace(".markDup.GCtagged", "")
pred_file_fragle_prl["FragleDetection"] = np.where(pred_file_fragle_prl["ctDNA_Burden"]>=0.01, 1, 0)
pred_file_fragle_nrh = pd.read_csv("/Users/alexandra/PhD/NeoRheaStudy/FragmentationPatterns/fragle/Fragle_Neorhea.csv", sep=",")
# pred_file_fragle_nrh = pd.concat([pred_file_fragle_nrh, fragle_valid], axis=0)
pred_file_fragle_nrh["FragleDetection"] = np.where(pred_file_fragle_nrh["ctDNA_Burden"]>=0.01, 1, 0)
pred_file_fragle_nrh["PatientId"] = pred_file_fragle_nrh["Sample_ID"].str.replace(".markDup.GCtagged", "")
pred_file_fragle_snr = pd.read_csv("/Users/alexandra/PhD/SynergyStudy/FragmentationPatterns/fragle/Fragle_Synergy.csv", sep=",")
# pred_file_fragle_snr = pd.concat([pred_file_fragle_snr, fragle_snr], axis=0)
pred_file_fragle_snr["PatientId"] = pred_file_fragle_snr["Sample_ID"].str.replace(".markDup.GCtagged", "")
pred_file_fragle_snr["FragleDetection"] = np.where(pred_file_fragle_snr["ctDNA_Burden"]>=0.01, 1, 0)

# test data
predictions_test_snr = pd.read_csv(pred_file_test, sep="\t")
predictions_test_snr["Predictedlabel"] = predictions_test_snr["Predictedlabel"]/100
predictions_test_snr = predictions_test_snr[((predictions_test_snr["PatientId"].str.contains("SNR")) | (predictions_test_snr["PatientId"].str.contains("NIPT-PIJB")))]
data_test_snr = pd.merge(predictions_test_snr, metadata, on="PatientId")
data_test_snr["ALFAssayDetection"] = np.where(data_test_snr["Predictedlabel"]>=0.03, 1, 0)
data_test_snr_all = pd.merge(data_test_snr, pred_file_fragle_snr,  on="PatientId")
data_test_snr_all["ALFAssayDetection"] = data_test_snr_all["ctDNADetected"]

# validation data
predictions_validation = pd.read_csv(pred_file_validation, sep="\t")
predictions_validation["Predictedlabel"] = predictions_validation["Predictedlabel"]/100
data_validation = pd.merge(predictions_validation, metadata, on="PatientId")
data_validation["ichorDetection"] = np.where(data_validation["ichorTF"]>=0.03, 1, 0)
data_validation["ALFAssayDetection"] = np.where(data_validation["Predictedlabel"]>=0.03, 1, 0)

predictions_validation_prl = data_validation[((data_validation["PatientId"].str.contains("PRL")) | (data_validation["PatientId"].str.contains("Genome-IJB")))]
data_validation_prl_all = pd.merge(predictions_validation_prl, pred_file_fragle_prl, on="PatientId")

predictions_validation_nrh = data_validation[((data_validation["PatientId"].str.contains("NRH")) | (data_validation["PatientId"].str.contains("Genome-IJB")))]
data_validation_nrh_all = pd.merge(predictions_validation_nrh, pred_file_fragle_nrh, on="PatientId")


## run classifyier

y_task = "Label"

test_snr=data_test_snr_all
# validation_prl=data_validation_prl_all
# validation_nrh=data_validation_nrh_all

validation_prl = data_validation_prl_all[data_validation_prl_all[y_task].notna()]
validation_nrh = data_validation_nrh_all[data_validation_nrh_all[y_task].notna()]


#ichor_features
x_feature_test=["ichorTF"]
y_feature_test=[y_task]
x_feature_validation=["ichorTF"]
y_feature_validation=[y_task]
probabilities_test_srn_ichor, probabilities_prl_ichor, probabilities_nrh_ichor, y_srn, y_prl, y_nrh = \
    perform_classification(test_snr, validation_prl, validation_nrh, x_feature_test, y_feature_test, x_feature_validation, y_feature_validation,"ichor")


#fragle_features
x_feature_test=["ctDNA_Burden"]
y_feature_test=[y_task]
x_feature_validation=["ctDNA_Burden"]
y_feature_validation=[y_task]
probabilities_test_srn_fragle, probabilities_prl_fragle, probabilities_nrh_fragle,  y_srn, y_prl, y_nrh = \
    perform_classification(test_snr, validation_prl, validation_nrh, x_feature_test, y_feature_test, x_feature_validation, y_feature_validation, "fragle")

#alfassay_features
x_feature_test=["Predictedlabel"]
y_feature_test=[y_task]
x_feature_validation=["Predictedlabel"]
y_feature_validation=[y_task]
probabilities_test_srn_alfassay, probabilities_prl_alfassay, probabilities_nrh_alfassay,  y_srn, y_prl, y_nrh = \
    perform_classification(test_snr, validation_prl, validation_nrh, x_feature_test, y_feature_test, x_feature_validation, y_feature_validation, "alfassay")

#alfassay_ichor_features
x_feature_test = ["Predictedlabel", "ichorTF"]
y_feature_test = [y_task]
x_feature_validation = ["Predictedlabel", "ichorTF"]
y_feature_validation = [y_task]
probabilities_test_srn_alfassay_ichor, probabilities_prl_alfassay_ichor, probabilities_nrh_alfassay_ichor,  y_srn, y_prl, y_nrh = \
    perform_classification(test_snr, validation_prl, validation_nrh, x_feature_test, y_feature_test, x_feature_validation, y_feature_validation, "alfassay_ichor")



# plot synergy
title="Classification models performance on Synergy healthyVScancer"
output_file = "/plots/classifier_t3_dagip/"
caption = '''Classification models performance on Synergy,\n ALFAssay model trained on SNR and  predicting TF.'''

V.plot_roc_fpr_tpr_multimodels(np.asarray(y_srn).squeeze(), np.asarray(probabilities_test_srn_ichor)[:,1],
                                 np.asarray(probabilities_test_srn_fragle)[:,1],
                                 np.asarray(probabilities_test_srn_alfassay)[:,1],
                                 np.asarray(probabilities_test_srn_alfassay_ichor)[:,1],
                                 output_file, title, 1, caption)

## plot pearl
title="Classification models performance on Pearl healthyVScancer"
output_file = "/plots/classifier_t3_dagip/"
caption = '''Classification models performance on Pearl data,\n ALFAssay model trained on SNR and  predicting TF.'''

V.plot_roc_fpr_tpr_multimodels(np.asarray(y_prl).squeeze(), np.asarray(probabilities_prl_ichor)[:,1],
                                 np.asarray(probabilities_prl_fragle)[:,1],
                                 np.asarray(probabilities_prl_alfassay)[:,1],
                                 np.asarray(probabilities_prl_alfassay_ichor)[:,1],
                                 output_file, title, 1, caption)


## plot neorhea
title="Classification models performance on Neorhea healthyVScancer" #-cancer/healthy
output_file = "/plots/classifier_t3_dagip/"
caption = '''Classification models performance on Neorhea data,\n ALFAssay model trained on SNR and  predicting TF.'''

V.plot_roc_fpr_tpr_multimodels(np.asarray(y_nrh).squeeze(), np.asarray(probabilities_nrh_ichor)[:,1],
                                 np.asarray(probabilities_nrh_fragle)[:,1],
                                 np.asarray(probabilities_nrh_alfassay)[:,1],
                                 np.asarray(probabilities_nrh_alfassay_ichor)[:,1],
                                 output_file, title, 1, caption)