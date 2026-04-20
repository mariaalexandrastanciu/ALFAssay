# Created by alexandra at 21/08/2023

from . import DataWrapperNN as dwnn
import time
from . import utils as u
import numpy as np
import torch as t
from . import Logger


class NNWrapper:
    def __init__(self, MainModel, losses):
        self.model = MainModel
        self.losses = losses
        self.dev = MainModel.dev

    def predict(self, input_data):
        self.model.eval()  # removes the noise inserted by dropouts

        # duplicates random input data such that batch size to be multiple of sample size
        input_data = dwnn.batch_padding(input_data, 1)
        loader = dwnn.dataWrapper(batch_size=1, num_workers=0, data=input_data, shuffle=False)

        pred_ys = []
        loss_dataset = []
        start_time = time.time()
        patientIds = []
        ctDNADetectionList = []
        accuracies = []
        bias_corrections = []
        corrected_xs = []
        true_xs = []
        # print("test param :", list(self.model.parameters()))
        for si, sample in enumerate(loader):
            fragmentSize, sampleNames, label, ctDNADetection, ctDNAbyVAF, ichorTF = sample
            task = ichorTF
            # correctFragmentSize = self.correct_data(fragmentSize)
            # ratio = self.calculate_ratio(correctFragmentSize)

            predicted,bias_correction, corrected_x = self.model(fragmentSize)
            # predicted = self.model(fragmentSize)
            truth = task.squeeze()

            # criterion = t.nn.MSELoss()
            # loss_tensor = criterion(predicted,  truth)
            criterion = u.RMSELoss
            loss_tensor = criterion(predicted, truth)
            #
            # loss = t.nn.MSELoss()
            # loss_tensor = loss(predicted, truth)


            pred_ys.append(predicted)  #[t.round(t.mean(predicted_arm_labels))]
            patientIds.append(sampleNames)
            ctDNADetectionList.append(truth)
            bias_corrections.append(bias_correction)
            corrected_xs.append(corrected_x)
            true_xs.append(fragmentSize.view(-1, 204))

            loss_dataset.append(float(loss_tensor.data.cpu().numpy()))
        accuracies = 0 #[u.accuracy_fn(ctDNADetectionList, pred_ys)]
        # print("\tprediction Time %s s" % round((time.time() - start_time), 2))
        return pred_ys, loss_dataset, patientIds, ctDNADetectionList,  accuracies , bias_corrections, corrected_xs, true_xs


    def fit(self, input_data, epochs=5, lr_scheduler_reduce=0.5, batch_size=1, save_model_every=10,
            weight_decay=1e-2, learning_rate=1e-3, LOG=True, learn_sum=False):

        if LOG:
            self.logger = Logger.MetaLogger(self.model, port=6001)

        ### data loader ###
        input_data = dwnn.batch_padding(input_data, batch_size)
        loader = dwnn.dataWrapper(batch_size=batch_size, num_workers=0, data=input_data)

        #######PRINT INFO##############
        parameters = list(self.model.parameters())
        used_params = []
        for i in parameters:
            if i.requires_grad:
                used_params += list(i.data.cpu().numpy().flat)
        # print('\tNumber of parameters=', len(used_params))

        ########OPTIMIZER##########
        optimizer = t.optim.Adam(self.model.parameters(), lr=learning_rate, weight_decay=weight_decay)
        scheduler = t.optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='min', factor=lr_scheduler_reduce,
                                                           patience=10000, verbose=True, threshold=1e-4,
                                                           threshold_mode='rel', cooldown=0, min_lr=1e-10, eps=1e-08)
        epoch_iteration = 0
        loss_value = 0.0
        t.autograd.set_detect_anomaly(True)

        all_losses = []
        all_pears_cor = []
        lrs = []
        last_epoch_loss = []
        while epoch_iteration < epochs:
            start = time.time()
            epoch_losses = []
            truth_list = []
            predicted_list =[]
            patientIds=[]
            # if epoch_iteration % 100 == 0:
            #     print("ok")

            last_epoch_loss = []
            for sample_iteration, sample in enumerate(loader):
                # see doc on true_y in dw.multiTaskCollate() function
                fragmentSize, sampleNames,  label, ctDNADetection, ctDNAbyVAF, ichorTF = sample
                task = ichorTF

                # correctFragmentSize = self.correct_data(fragmentSize)
                # ratio = self.calculate_ratio(correctFragmentSize)

                predicted, bias_correction, corrected_x = self.model(fragmentSize)   # calls forward
                # predicted, bias_correction, corrected_x = self.model(fragmentSize)   # calls forward

                truth = task.squeeze()
                # criterion = t.nn.MSELoss()
                # loss_tensor = criterion(predicted,  truth)
                criterion = u.RMSELoss
                loss_tensor = criterion(predicted, truth)

                # loss_tensor = criterion(predicted,  truth)
                optimizer.zero_grad()
                epoch_losses += [float(loss_tensor.data.cpu().numpy())]

                loss_tensor.backward()

                last_epoch_loss.append(float(loss_tensor.data.cpu().numpy()))
                patientIds.append(sampleNames)

                if LOG:
                    self.logger.update_weights(epoch_iteration)

                optimizer.step()
                optimizer.zero_grad()

                truth_list.append(truth)
                # predicted_original = self.target_transform.inverse_transform(predicted)

                predicted_list.append(predicted)

            if epoch_iteration % 100 == 0:
                # print(np.mean(epoch_losses))
                print(" train epoch ", epoch_iteration, "; loss mean=", np.mean(epoch_losses))
                # print("\ntruth train: \n" , truth_list)
                # print("\npredicted train: \n",  predicted_list)
                # print(" epoch ", epoch_iteration, "time", round(end - start, 2))
                # print("train param :", list(self.model.parameters()))
            end = time.time()
            scheduler.step(np.mean(epoch_losses))
            if LOG:
                self.logger.writeTensorboardLog(epoch_iteration, np.mean(epoch_losses), end - start,
                                                self.model.embeds(t.arange(0, 20)))

            epoch_iteration += 1
            all_losses += [np.mean(epoch_losses)]
        # print(" last train epoch ", epoch_iteration, "; loss mean=", np.mean(last_epoch_loss))
        # print("last losses on train", last_epoch_loss)
        # print("\ntruth train: \n", truth_list)
        # print("\npredicted train: \n", predicted_list)
        # print(" epoch ", epoch_iteration, "time", round(end - start, 2))
        # print("train param :", list(self.model.parameters()))
        # print("patients", patientIds)

    def get_params(self, deep=False):
        return {}

    def get_model(self):
        return self.model
