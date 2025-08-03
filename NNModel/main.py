import pandas as pd
import preprocessing_data as ppd

import os
import torch as t
# os.environ['KMP_DUPLICATE_LIB_OK'] = 'True'
import numpy as np
import NNModel
from NNWrapper import NNWrapper
from sklearn.model_selection import train_test_split
from sklearn.model_selection import KFold, StratifiedKFold
import time
import ApplicationSetup as AS
import CONSTANTS as c
# import DataWrapperNN as DWNN
import CustomLosses as losses
import ModelParameters as MP
import Visualisation as V
dev = "cpu"  # AS.if_cuda(); TODO: should be fixed
AS.application_setup()

import pickle

tasks = [c.Label, c.ctDNADetected, c.VAFg0p001]
task = c.ctDNADetected
try_model = "_t36"
def main():
    # input_file = "/Users/alexandra/PhD/FragmentationPatterns/Data/NN_5MB/NN_5MB_scaled.csv"
    # test_studies = ["PearlStudy",  "healthy_sWGS"]
    # validation_studies = ["NeoRheaStudy", "healthy_sWGS"]
    # test_data, validation_data, test_labels, validation_labels,  fragmentSizeFeatureTest, fragmentSizeFeatureValidation \
    #     = ppd.get_data_for_wrapper(input_file, test_studies, validation_studies, task)

    # input_file = "/Users/alexandra/PhD/FragmentationPatterns/Data/NN_5MB/NN_5MB_pivot.csv"
    input_file = "/Users/alexandra/PhD/FragmentationPatterns/Data/NN_5MB/NN_5MB_pivot_with_ratios_.csv"
    test_studies = ["SNR", "NIPT-PIJB", "NRH"] #"NIPT-PIJB"
    validation_studies = ["PRL", "Genome-IJB"]
    # features = ["noReadsGCParagonFrag", "shortOverLongGCPara"]
    features = ["noReadsGCParagonFrag"]

    #train on different cohorts and validate on different cohorts
    wrapperDataTraining, wrapperDataValidation, labels_training, features_training, labels_validation, \
    features_validation, pt = ppd.get_data_wrapper_2(input_file, "ichorTF")

    inSilico_input_file = "/Users/alexandra/PhD/PyCharmProjects/ALFAssay/validations/data/synthetic_ctDNA_combined_all_fractions.csv"
    wrapperDataValidationInSilico,  labels_validationInSilico, features_validationInSilico = \
        ppd.get_data_insilico(inSilico_input_file, "ichorTF")

    # wrapperDataTraining, wrapperDataValidation, labels_training, features_training, labels_validation, \
    # features_validation, pt = ppd.get_nn_bias_corrected_data("ichorTF", features)

    #mixed metastatic
    # wrapperDataTraining, wrapperDataValidation, labels_training, features_training, labels_validation, \
    # features_validation = ppd.get_data_wrapper_mixed(input_file,  "ctDNADetected", features)


    # train on SRN positive and healthy
    # validation_studies = ["Neorhea", "Pearl"]
    # wrapperDataTraining, wrapperDataValidation, labels_training, features_training, labels_validation, \
    # features_validation = ppd.get_data_wrapper_synergy_train(input_file, "ctDNADetected", features, validation_studies)


    ## use dagip corrected fragements
    # dagip_file = "/Users/alexandra/PhD/PyCharmProjects/DAGIP/data/corrected_data_dagip.csv"
    # wrapperDataTraining, wrapperDataValidation, labels_training, features_training, labels_validation, \
    # features_validation, pt = ppd.get_data_wrapper_corrected(dagip_file, test_studies, validation_studies, "ichorTF")

    losses_global = [losses.MSE()]
    bins = 204 # 2141 # 204
    MainModel = NNModel.NeuralNetworkModel(input_dim=1*bins, output_dim=1,
                                           hidden_layers=40, dev=dev, batch_size=MP.batch_size)

    # dataset = np.array(data)
    kf = KFold(n_splits=MP.cv, shuffle=True)   # rmsd = []    mse = []
    nr_epochs = np.arange(800, 4001, 800)
    # nr_epochs = np.arange(200, 1401, 200)
    # nr_epochs = np.arange(1, 3, 1)
    # x_train, x_test = train_test_split(wrapperDataTraining, test_size=0.3, random_state=123, stratify=labels_training[:, task])

    # X = scaler.fit_transform(X)
    # X_v = scaler.transform(X_v)

    wrapper = NNWrapper(MainModel, losses_global, pt)  # builds the wrapper of the NN
    # mean_mse_train = []
    # mean_mse_test = []
    # mse_train = []
    # mse_test = []
    # true_ys_test =[]
    # preds_y_test=[]
    # x_train = [x_train[0], x_train[1], x_train[2], x_train[3], x_train[4], x_train[5], x_train[6], x_train[7], x_train[8], x_train[9],
    #            x_train[10], x_train[11], x_train[12], x_train[13], x_train[14], x_train[15], x_train[16], x_train[17],  x_train[18], x_train[19],
    #            x_train[20],  x_train[21], x_train[22], x_train[23],  x_train[24], x_train[25],
    #            x_train[26],  x_train[27], x_train[28], x_train[29],  x_train[30], x_train[31],
    #            x_train[32],  x_train[33], x_train[34], x_train[35],  x_train[36], x_train[37]]
    #
    # x_train = [x_train[0], x_train[1], x_train[6], x_train[7], x_train[13], x_train[14], x_train[15]]
               # x_train[16], x_train[17], x_train[18], x_train[19]]
    # x_train = [x_train[0],x_train[1], x_train[2] , x_train[18] x_train[19]]
    # x_train = [x_train[0]]
    #
    # wrapperDataTraining = [wrapperDataTraining[0], wrapperDataTraining[1], wrapperDataTraining[2], wrapperDataTraining[3], wrapperDataTraining[4],
    #                        wrapperDataTraining[5], wrapperDataTraining[6], wrapperDataTraining[7], wrapperDataTraining[8], wrapperDataTraining[9]]


    loss_on_train = []
    loss_on_test = []
    for epoch in nr_epochs:
        cv_start = time.time()
        mean_mse_train = []
        mean_mse_test = []
        mse_train = []
        mse_test = []
        true_ys_test = []
        preds_y_test = []
        patientIds_test = []
        bias_corrections = []
        corrected_xs = []
        true_xs = []
        for train, test in kf.split(wrapperDataTraining):

            start = time.time()
            trainData = [wrapperDataTraining[data_train] for data_train in train]    # training structures
            testData = [wrapperDataTraining[data_test] for data_test in test]
            wrapper.fit(trainData, epochs=epoch, lr_scheduler_reduce=MP.lr_scheduler_reduce,
                        batch_size=MP.batch_size, save_model_every=MP.save_model_sec, weight_decay=MP.weight_decay,
                        learning_rate=MP.learning_rate, LOG=MP.LOG, learn_sum=MP.learn_sum)
            end = time.time()
            # print("### time to train: ", end - start, " for ", epoch, "epochs")
            # prediction on train dataset

            ypred_train, loss_dataset_train, patientId_train, ctDNADetection_train,  accuracy_train, bias_correction_train, corrected_x_train, true_xs_train \
                = wrapper.predict(trainData, batch_size=MP.batch_size)

            # ypred_train, loss_dataset_train, patientId_train, ctDNADetection_train, accuracy_train\
            #     = wrapper.predict(trainData, batch_size=MP.batch_size)

            # our actual prediction
            # start = time.time()
            # pred_y_test, loss_dataset_test, patientId_test, true_y_test, accuracy_test, bias_correction, corrected_x, true_x = \
            #     wrapper.predict(testData, batch_size=MP.batch_size)

            pred_y_test, loss_dataset_test, patientId_test, true_y_test, accuracy_test, bias_correction, corrected_x, true_x = \
                wrapper.predict(testData, batch_size=MP.batch_size)
            # end = time.time()
            # print("### time to predict: ", end - start, " for ", epoch, "epochs")
            mse_train.append(np.mean(loss_dataset_train))
            mse_test.append(np.mean(loss_dataset_test))
            patientIds_test += patientId_test
            preds_y_test += pred_y_test
            true_ys_test += true_y_test
            bias_corrections += bias_correction
            corrected_xs += corrected_x
            true_xs += true_x
            mean_mse_train.append(np.mean(mse_train))
            mean_mse_test.append(np.mean(mse_test))

        cv_end = time.time()
        loss_on_train.append(np.mean(mean_mse_train))
        loss_on_test.append(np.mean(mean_mse_test))

        print("### epoch ", epoch, "\nLoss on Train:", np.mean(mean_mse_train), "\nLoss on Test",  np.mean(mean_mse_test))
        print("## time for CV:" , cv_end - cv_start, " for ", epoch, "epochs")
        print("\npatientId_test\n", patientIds_test)
        print("\nypred_test\n", preds_y_test)
        print("\nctDNADetection_test\n", true_ys_test)


    print("Mean loss function on Train:", np.mean(loss_on_train), "\n Mean loss on Test", np.mean(loss_on_test))
    title = "Loss function value over the epochs on training and on test data sets " + try_model
    V.draw_result(nr_epochs, loss_on_train, loss_on_test, title)
    title="train SNR and healthy ichorTF "
    output_file = "plots/"
    with t.no_grad():
        patientIds_flat = sum(sum(patientIds_test, []), [])
        df_results = pd.DataFrame(np.asarray([np.asarray(patientIds_flat), np.asarray(preds_y_test),
                                              np.asarray(true_ys_test)]).T,
                                  columns=["PatientId", "Predictedlabel", "Truelabel"])
        df_results.to_csv("results/results_train_ichorTF"+try_model + ".csv", sep="\t", index=False)

        df_correction_test = pd.DataFrame(np.asarray(t.stack(corrected_xs)).squeeze().T,
                                     columns=np.asarray(patientIds_flat))
        df_correction_test.to_csv("results/results_test_correction"+try_model + ".csv", sep="\t", index=False)

        df_bias_test = pd.DataFrame(np.asarray(t.stack(bias_corrections)).squeeze().T,
                               columns=np.asarray(patientIds_flat))
        df_bias_test.to_csv("results/results_test_bias"+try_model + ".csv", sep="\t", index=False)

        df_real_test = pd.DataFrame(np.asarray(t.stack(true_xs)).squeeze().T,
                                    columns=np.asarray(patientIds_flat))
        df_real_test.to_csv("results/results_real_bias"+try_model + ".csv", sep="\t", index=False)

        # V.plot_roc_fpr_tpr(ctDNADetection_test, ypred_test, output_file, title, 1)

    #save model
    with open('/Users/alexandra/PhD/PyCharmProjects/ALFAssayNN/models/ALFAssay' +try_model + ".pkl", 'wb') as f:
        pickle.dump(wrapper.get_model(), f)

   # validation
    pred_y_validation, loss_dataset_validation, patientId_validation, true_y_validation, accuracy_validation, bias_corrections_v, corrected_xs_v, true_xs_v = \
        wrapper.predict(wrapperDataValidation, batch_size=MP.batch_size)

    print( "\nLoss on validation:", np.mean(loss_dataset_validation))
    print("\nypred_validation\n", pred_y_validation)
    print("\nctDNADetection_validation\n", true_y_validation)

    # validation_in_silico
    pred_y_validation_in_silico, loss_dataset_validation_in_silico, patientId_validation_in_silico, \
    true_y_validation_in_silico, accuracy_validation_in_silico, bias_corrections_v_in_silico, corrected_xs_v_in_silico, \
    true_xs_v_in_silico = \
        wrapper.predict(wrapperDataValidationInSilico, batch_size=MP.batch_size)

    print("\nLoss on validation in silico:", np.mean(loss_dataset_validation_in_silico))
    # print("\nypred_validation\n", pred_y_validation)
    # print("\nctDNADetection_validation\n", true_y_validation)

    title = "validation NRH, PRL,Healthy ichorTF"
    output_file = "plots/"
    with t.no_grad():
        # fpr, tpr, _ = roc_curve(np.asarray(ctDNADetection_test), np.asarray(ypred_test))
        #
        patientIds_validation = sum(sum(patientId_validation, []), [])
        df_results = pd.DataFrame(np.asarray([np.asarray(patientIds_validation), np.asarray(pred_y_validation),
                                              np.asarray(true_y_validation)]).T,
                                  columns=["PatientId", "Predictedlabel", "Truelabel"])
        df_results.to_csv("results/results_validation_ichorTF"+try_model + ".csv", sep="\t", index=False)

        ### in silico
        patientIds_validation_in_silico = sum(sum(patientId_validation_in_silico, []), [])
        df_results = pd.DataFrame(np.asarray([np.asarray(patientIds_validation_in_silico), np.asarray(pred_y_validation_in_silico),
                                              np.asarray(true_y_validation_in_silico)]).T,
                                  columns=["PatientId", "Predictedlabel", "Truelabel"])
        df_results.to_csv("results/results_validationInSilico_ichorTF" + try_model + ".csv", sep="\t", index=False)


        df_correction = pd.DataFrame(np.asarray(t.stack(corrected_xs_v)).squeeze().T,
                                  columns=np.asarray(patientIds_validation))
        df_correction.to_csv("results/results_validation_correction_"+try_model + ".csv", sep="\t", index=False)

        df_bias = pd.DataFrame(np.asarray(t.stack(bias_corrections_v)).squeeze().T,
                                     columns=np.asarray(patientIds_validation))
        df_bias.to_csv("results/results_validation_bias"+try_model + ".csv", sep="\t", index=False)

        df_real = pd.DataFrame(np.asarray(t.stack(true_xs_v)).squeeze().T,
                               columns=np.asarray(patientIds_validation))
        df_real.to_csv("results/results_validation_real"+try_model + ".csv", sep="\t", index=False)

        # # Plotting
        # plt.plot(fpr, tpr, marker='.')
        # plt.ylabel('True Positive Rate')
        # plt.xlabel('False Positive Rate')
        # plt.show()
        # V.plot_roc_fpr_tpr(ctDNADetection_validation, ypred_validation, output_file, title, 1)

    print('The monkeys are listening')


if __name__ == '__main__':
    import sys

    sys.exit(main())
