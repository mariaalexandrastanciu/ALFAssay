# Created by alexandra at 31/10/2025
import pandas as pd
from . import preprocessing_data as ppd

import torch as t
import numpy as np
from . import NNModel
from . import NNWrapper
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.metrics import mean_squared_error, r2_score
from scipy.stats import spearmanr
import time
from . import ApplicationSetup as AS
from . import CONSTANTS as c
from . import CustomLosses as losses
from . import ModelParameters as MP
from . import Visualisation as V
import pickle

dev = "cpu"  # AS.if_cuda();
AS.application_setup()

tasks = [c.Label, c.ctDNADetected, c.VAFg0p001]
task = c.ctDNADetected  # <-- used only for labels indexing; regression target is ichorTF below
try_model = "_tPLOS1"

def to_numpy_flat(x):
    import numpy as np
    import torch as t

    def _as_1d_np(v):
        # torch tensor
        if isinstance(v, t.Tensor):
            return v.detach().cpu().numpy().reshape(-1)
        # numpy array/scalar
        if isinstance(v, np.ndarray):
            return v.reshape(-1)
        # python scalar
        if np.isscalar(v):
            return np.array([v], dtype=float)
        # something else iterable -> recurse
        if isinstance(v, (list, tuple)):
            parts = []
            for vi in v:
                arr = _as_1d_np(vi)
                if arr.size:
                    parts.append(arr)
            if len(parts) == 0:
                return np.array([], dtype=float)
            return np.concatenate(parts, axis=0)
        # fallback
        return np.array([v], dtype=float)

    out = _as_1d_np(x)
    return out.astype(float, copy=False)  # ensure numeric dtype

def main():
    # ---------------------------
    # Load training/validation data
    # ---------------------------
    input_file = "FeatureEngineering/data/NN_5MB_features.csv"


    # returns:
    # wrapperDataTraining, wrapperDataValidation, labels_training, features_training, labels_validation, features_validation
    wrapperDataTraining, wrapperDataValidation, labels_training, features_training, labels_validation, \
        features_validation = ppd.get_data_wrapper(input_file, "ichorTF")

    # ---------------------------
    # Prepare arrays for grouped + stratified CV
    # ---------------------------
    # labels_training columns are: ["PatientId", "Label", "ctDNADetected", "VAF", "ichorTF", "study"]
    patient_ids = labels_training[:, 0].astype(str)
    studies = labels_training[:, 5].astype(str)

    # Regression target: ichorTF (continuous)
    y_cont = labels_training[:, 4].astype(float)

    assert len(wrapperDataTraining) == len(patient_ids) == len(studies) == len(y_cont), \
        "Length mismatch between data and metadata."

    # Build binned target for stratification. We use global quantile bins but
    # cap the number of bins by the smallest study size to avoid empty bins per fold.
    from collections import Counter
    n_splits = MP.cv
    study_counts = Counter(studies)
    min_per_study = min(study_counts.values())
    max_bins_by_counts = max(1, min(10, min_per_study // n_splits))
    if max_bins_by_counts > 1:

        y_bins = pd.qcut(y_cont, q=max_bins_by_counts, labels=False, duplicates='drop')
    else:
        y_bins = np.zeros_like(y_cont, dtype=int)

    # Composite stratification label: study × y-bin
    strat_labels = np.array([f"{s}_{b}" for s, b in zip(studies, y_bins)])

    # ---------------------------
    # Cross-validation setup (precompute splits once)
    # ---------------------------
    cv = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=42)
    splits = list(cv.split(X=np.zeros_like(y_bins), y=strat_labels, groups=patient_ids))

    # ---------------------------
    # NN / training setup
    # ---------------------------
    losses_global = [losses.MSE()]
    bins = 408  # input_dim = bins * channels
    nr_epochs = np.arange(1200, 4801, 800)
    # nr_epochs = np.arange(800, 801, 800)

    # ---------------------------
    # CV over epochs (no leakage: new model + new PT per fold)
    # ---------------------------
    loss_on_train = []
    loss_on_test = []

    for epoch in nr_epochs:
        cv_start = time.time()

        fold_mse_train = []
        fold_mse_test = []

        # Per-epoch collectors (OPTIONAL detailed logging)
        all_patientIds_test = []
        all_preds_y_test = []
        all_true_y_test = []
        all_bias_corrections = []
        all_corrected_xs = []
        all_true_xs = []

        for fold_idx, (train_idx, test_idx) in enumerate(splits, start=1):
            # Sanity checks
            assert len(set(studies[test_idx])) >= 3, f"[Fold {fold_idx}] Missing a study in test: {set(studies[test_idx])}"
            assert set(patient_ids[train_idx]).isdisjoint(set(patient_ids[test_idx])), "[Fold {fold_idx}] Patient leakage detected."

            # Build fold data
            trainData = [wrapperDataTraining[i] for i in train_idx]
            testData  = [wrapperDataTraining[i] for i in test_idx]


            # Fresh model & wrapper per fold (avoid cross-fold training carryover)
            MainModel = NNModel.NeuralNetworkModel(
                input_dim=1 * bins, output_dim=1,
                hidden_layers=40, dev=dev, batch_size=MP.batch_size
            )
            wrapper = NNWrapper.NNWrapper(MainModel, losses_global)

            # Train
            wrapper.fit(
                trainData, epochs=epoch, lr_scheduler_reduce=MP.lr_scheduler_reduce,
                batch_size=MP.batch_size, save_model_every=MP.save_model_sec,
                weight_decay=MP.weight_decay, learning_rate=MP.learning_rate,
                LOG=MP.LOG, learn_sum=MP.learn_sum
            )

            # Predict on train/test
            ypred_train, loss_dataset_train, patientId_train, ytrue_train, accuracy_train, bias_correction_train, corrected_x_train, true_x_train = \
                wrapper.predict(trainData)

            ypred_test, loss_dataset_test, patientId_test, ytrue_test, accuracy_test, bias_correction_test, corrected_x_test, true_x_test = \
                wrapper.predict(testData)

            # Regression metrics (fold-level)
            # Note: your wrapper returns per-sample losses; we also compute RMSE/R2 for monitoring
            ytrue_test_np = to_numpy_flat(ytrue_test).ravel()
            ypred_test_np = to_numpy_flat(ypred_test).ravel()

            rmse_te = np.sqrt(mean_squared_error(ytrue_test_np, ypred_test_np))
            r2_te = r2_score(ytrue_test_np, ypred_test_np) if len(np.unique(ytrue_test_np)) > 1 else np.nan
            rho_te, _ = spearmanr(ytrue_test_np, ypred_test_np) if len(ytrue_test_np) > 1 else (np.nan, None)

            # Aggregate losses (your original code uses mean MSE per fold)
            fold_mse_train.append(np.mean(loss_dataset_train))
            fold_mse_test.append(np.mean(loss_dataset_test))

            # Collect for optional export per epoch
            all_patientIds_test += patientId_test
            all_preds_y_test += ypred_test
            all_true_y_test += ytrue_test
            all_bias_corrections += bias_correction_test
            all_corrected_xs += corrected_x_test
            all_true_xs += true_x_test

            print(f"[Epoch {epoch} | Fold {fold_idx}/{n_splits}] "
                  f"MSE(train)={np.mean(loss_dataset_train):.6f}  "
                  f"MSE(test)={np.mean(loss_dataset_test):.6f}  "
                  f"RMSE(test)={rmse_te:.6f}  R2(test)={r2_te:.4f}  Spearman={rho_te:.4f}")

        cv_end = time.time()

        mean_mse_train = np.mean(fold_mse_train)
        mean_mse_test = np.mean(fold_mse_test)
        loss_on_train.append(mean_mse_train)
        loss_on_test.append(mean_mse_test)

        print("### epoch ", epoch,
              "\nLoss on Train:", mean_mse_train,
              "\nLoss on Test",  mean_mse_test)
        print("## time for CV:", cv_end - cv_start, " for ", epoch, "epochs")

        # Flatten nested patient IDs (wrapper returns nested lists)
        with t.no_grad():
            try:
                patientIds_flat = sum(sum(all_patientIds_test, []), [])
            except Exception:
                # If already flat, keep as-is
                patientIds_flat = all_patientIds_test

            df_results = pd.DataFrame(
                np.asarray([np.asarray(patientIds_flat), np.asarray(all_preds_y_test),
                            np.asarray(all_true_y_test)]).T,
                columns=["PatientId", "PredictedValue", "TrueValue"]
            )
            df_results.to_csv(f"NNModel/results/results_cv_epoch_{epoch}{try_model}.csv", sep="\t", index=False)

    # ---------------------------
    # Plot CV loss over epochs
    # ---------------------------
    print("Mean loss function on Train:", np.mean(loss_on_train),
          "\n Mean loss on Test", np.mean(loss_on_test))
    title = "Loss over epochs on training vs test (CV, grouped+stratified) " + try_model
    V.draw_result(nr_epochs, loss_on_train, loss_on_test, title)

    # ---------------------------
    # Train a FINAL model on ALL training data (for export + external validation)
    # ---------------------------
    best_epoch = int(nr_epochs[np.argmin(loss_on_test)])  # pick epoch with lowest mean CV test loss
    print(f"Training final model on full training set for {best_epoch} epochs...")


    MainModel_final = NNModel.NeuralNetworkModel(
        input_dim=1 * bins, output_dim=1,
        hidden_layers=40, dev=dev, batch_size=MP.batch_size
    )
    wrapper_final = NNWrapper.NNWrapper(MainModel_final, losses_global)

    wrapper_final.fit(
        wrapperDataTraining, epochs=best_epoch, lr_scheduler_reduce=MP.lr_scheduler_reduce,
        batch_size=MP.batch_size, save_model_every=MP.save_model_sec,
        weight_decay=MP.weight_decay, learning_rate=MP.learning_rate,
        LOG=MP.LOG, learn_sum=MP.learn_sum
    )

    # Save the final model
    with open('NNModel/models/ALFAssay' + try_model + ".pkl", 'wb') as f:
        pickle.dump(wrapper_final.get_model(), f)

    # ---------------------------
    # External validation on held-out cohorts
    # ---------------------------
    pred_y_validation, loss_dataset_validation, patientId_validation, true_y_validation, accuracy_validation, \
        bias_corrections_v, corrected_xs_v, true_xs_v = \
        wrapper_final.predict(wrapperDataValidation)

    print("\nLoss on validation:", np.mean(loss_dataset_validation))
    # print("\nypred_validation\n", pred_y_validation)
    # print("\ntrue_validation\n", true_y_validation)


    # Save validation outputs
    with t.no_grad():
        try:
            patientIds_validation = sum(sum(patientId_validation, []), [])
        except Exception:
            patientIds_validation = patientId_validation

        df_results = pd.DataFrame(np.asarray([np.asarray(patientIds_validation), np.asarray(pred_y_validation),
                                              np.asarray(true_y_validation)]).T,
                                  columns=["PatientId", "PredictedValue", "TrueValue"])
        df_results.to_csv("NNModel/results/results_validation_ichorTF" + try_model + ".csv", sep="\t", index=False)

    print('The monkeys are listening')
    return 0


if __name__ == '__main__':
    import sys
    sys.exit(main())
