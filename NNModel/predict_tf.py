# Created by alexandra at 07/04/2026
#!/usr/bin/env python3
# Created by alexandra
import argparse
import pickle
import pandas as pd

from . import preprocessing_data as ppd
from . import NNWrapper
from . import CustomLosses as losses
import torch as t


def to_flat_list(x):
    try:
        return sum(sum(x, []), [])
    except Exception:
        return x


def main():
    parser = argparse.ArgumentParser(description="Predict tumor fraction from a feature file.")
    parser.add_argument("input_file", help="Path to the feature file")
    parser.add_argument(
        "--model",
        default="NNModel/models/ALFAssay_model.pkl",
        help="Path to the saved model bundle"
    )
    parser.add_argument(
        "--output",
        default="NNModel/results/predictions.tsv",
        help="Path to the output prediction file"
    )
    args = parser.parse_args()

    # Prepare data using the same preprocessing function as training
    # Assumes the file can be passed through the same wrapper preparation logic
    wrapperDataTraining, wrapperDataValidation, labels_training, features_training, labels_validation, features_validation = \
        ppd.get_data_wrapper(args.input_file, "ichorTF")

    # For prediction-only usage, we take the prepared dataset that corresponds
    # to the provided input file. If your helper returns everything in the
    # training slot for standalone files, switch to wrapperDataTraining instead.
    data_to_predict = wrapperDataValidation
    labels_to_use = labels_validation

    # Load model bundle
    with open(args.model, "rb") as f:
        model = pickle.load(f)

    # Expected saved object format:
    # {"model": trained_model, "power_transformer": fitted_power_transformer}

    wrapper = NNWrapper.NNWrapper(model, [losses.MSE()])

    pred_y, loss_dataset, patient_ids, true_y, _, _, _, _ = wrapper.predict(data_to_predict)

    patient_ids = to_flat_list(patient_ids)
    with t.no_grad():
        df = pd.DataFrame({
            "PatientId": patient_ids,
            "PredictedValue": pred_y,
            "TrueValue": true_y
        })

        df.to_csv(args.output, sep="\t", index=False)
    print(f"Predictions saved to: {args.output}")


if __name__ == "__main__":
    raise SystemExit(main())