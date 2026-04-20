# Created by alexandra at 18/12/2023
import numpy as np
import pandas as pd



def get_data_wrapper(file, task):
    metadata = pd.read_csv("refdata/AllMetaData.csv", sep="\t")
    labels_cols = ["PatientId", "Label", "ctDNADetected", "VAF", "ichorTF", "study"]
    data = pd.read_csv(file, sep="\t")
    data = data.fillna(0)

    data_feat = data
    metadata = metadata[metadata[task].notna()]
    metadata["ctDNADetected"] = metadata["ctDNADetected"].apply(float)
    metadata["ichorTF"] = metadata["ichorTF"].apply(float)
    metadata["ichorTF"] = np.where(metadata["ichorTF"]<0.03, 0 ,metadata["ichorTF"]*100)
    metadata["Label"] = metadata["Label"].apply(float)
    metadata["VAF"] = metadata["VAF"].apply(float)

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
        totalReads = [patient_data[patient_data["Feature"] == "noReadsGCParagonFrag"][patient].apply(pd.to_numeric).values]
        labels = metadata_patient[labels_cols].values
        wrapperDataTraining.append([np.asarray([shortRatio, totalReads]).T, labels])
        # wrapperDataTraining.append([np.asarray(shortRatio).T, labels])
        labels_training.append(labels)
        features_training.append(np.asarray([shortRatio, totalReads]).T)
        # features_training.append(np.asarray(shortRatio).T)

    wrapperDataValidation = []
    labels_validation = []
    features_validation = []
    for patient in validation_columns:
        patient_data = data_feat[["Feature", patient]]
        metadata_patient = metadata[metadata["PatientId"] == patient]
        shortRatio = [patient_data[patient_data["Feature"] == "ratioShort"][patient].apply(pd.to_numeric).values]
        longRatio = [patient_data[patient_data["Feature"] == "ratioLong"][patient].apply(pd.to_numeric).values]
        totalReads = [
            patient_data[patient_data["Feature"] == "noReadsGCParagonFrag"][patient].apply(pd.to_numeric).values]
        labels = metadata_patient[labels_cols].values
        wrapperDataValidation.append([np.asarray([shortRatio, totalReads]).T, labels])
        # wrapperDataValidation.append([np.asarray(shortRatio).T, labels])
        labels_validation.append(labels)
        features_validation.append(np.asarray([shortRatio, totalReads]).T)
        # features_validation.append(np.asarray(shortRatio).T)

    labels_training = np.asarray(labels_training).squeeze()
    features_training = np.asarray(features_training).squeeze()
    labels_validation = np.asarray(labels_validation).squeeze()
    features_validation = np.asarray(features_validation).squeeze()


    return wrapperDataTraining, wrapperDataValidation, labels_training, features_training, labels_validation,\
           features_validation
