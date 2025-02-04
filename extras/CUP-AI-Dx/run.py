#!/usr/bin/env ipython
import argparse
import json
import os
import time as tt
from pathlib import Path

import numpy as np
import pandas as pd
from keras.models import load_model
from sklearn.metrics import (
    classification_report,
    confusion_matrix,
    f1_score,
    precision_score,
    recall_score,
)

from datasets import load_dataset

# Model expects a matrix of gene x samples where genes are in enterez ids and
# the values are log2(TPM + 1) data
from utils import reshape_data_1d


def load_label_encodings(data_dir):
    with open(f"{data_dir}/label_encoding.json") as fin:
        _encoding = json.load(fin)
    int_to_label = dict((int(k), v) for k, v in _encoding.items())
    label_to_int = dict((v, int(k)) for k, v in _encoding.items())
    return int_to_label, label_to_int


def sk_report_performance(labels, predictions, prefix: str = "average_"):
    metrics = ["precision", "recall", "f1_score"]
    metrics = list(map(lambda x: f"{prefix}{x}", metrics))
    scores = [
        precision_score(labels, predictions, average="weighted", zero_division=0),
        recall_score(labels, predictions, average="weighted", zero_division=0),
        f1_score(labels, predictions, average="weighted"),
    ]
    df = pd.DataFrame({"metric": metrics, "score": scores})
    return df


def predict(
    input: pd.DataFrame, labels: pd.Series, model, verbose: bool = False
) -> dict:
    # input is samples x features
    samples = input.index
    reshaped_x = reshape_data_1d(input.values)  # add a new dimension
    if verbose:
        print("INPUT")
        print(input)
        print("RESHAPED")
        print(reshaped_x.shape)
    pred_prob = model.predict(reshaped_x)  # samples x features
    pred_labels = np.argmax(pred_prob, axis=1)  # labels for samples

    pred_labels = [from_int_to_label[k] for k in pred_labels]
    with_pred = pd.DataFrame({"prediction": pred_labels, "sample": samples})

    result: dict = {
        "cm": None,
        "report": None,
        "prediction": with_pred,
        "metrics": None,
    }
    if labels.any():
        cm_labels = list(set(pred_labels) | set(labels))
        cm_labels.sort()
        _confusion_matrix = confusion_matrix(labels, pred_labels, labels=cm_labels)
        cm = pd.DataFrame(_confusion_matrix, index=cm_labels, columns=cm_labels)
        report = classification_report(labels, pred_labels, labels=cm_labels)
        result["report"] = report
        result["cm"] = cm
        metrics = sk_report_performance(labels, pred_labels)
        print(metrics)
        result["metrics"] = metrics
        acc = np.sum(_confusion_matrix.diagonal()) / np.sum(_confusion_matrix)
        print(f"overall accuracy: {acc * 100}%")
    return result


def load_input(file: str, data_dir: str, label_column: str):
    """`file` is either a gene x sample matrix, or a sample x gene matrix
        the latter case is if labels are known (in a column called "tumor_type")
        and the file is provided for testing purposes.

    entries are gene expression values in log2 tpm

    :returns: a sample x gene matrix (which the model expects),
        filtered for cup-ai-dx's selected genes
        and a series of labels if provided
    """
    df = pd.read_csv(file, index_col=0)
    if label_column in df.columns:
        labels: pd.Series = df[label_column]
        df = np.transpose(df.drop([label_column], axis=1))
    else:
        labels = pd.Series(None, index=df.columns, name=label_column)

    target_genes = pd.read_csv(
        f"{data_dir}/features_791.csv", header=None, index_col=0
    ).index

    if "X" not in str(df.index[0]):
        df.index = [f"X{i}" for i in df.index]

    rows = []
    for t in target_genes:
        try:
            rows.append(df.loc[[t],])
        except KeyError:
            rows.append(pd.DataFrame(0, columns=df.columns, index=[t]))
    selected = pd.concat(rows, axis=0)
    transposed: pd.DataFrame = np.transpose(selected)
    return transposed, labels


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--data_dir", default="data", action="store")
    parser.add_argument(
        "-i",
        "--input",
        default="",
        required=True,
        help="csv file containing expression data",
        action="store",
    )
    parser.add_argument("-p", "--prefix", default="cupai", help="Prefix for output")
    parser.add_argument("-o", "--outdir", default=".", help="Output directory")
    parser.add_argument(
        "-l",
        "--label_column",
        default="label",
        help="For test datasets (sample x genes), the column indicating the label of the sample",
        action="store",
    )
    parser.add_argument(
        "--model",
        required=False,
        default="inception",
        choices=["inception", "cnn", "resnet"],
    )
    parser.add_argument(
        "--models_dir", default="models", action="store", required=False
    )
    args = parser.parse_args()
    from_int_to_label, from_label_to_int = load_label_encodings(args.data_dir)

    start_time = tt.time()

    model_files = {
        "inception": f"{args.models_dir}/inception_net_1d.h5",
        "cnn": f"{args.models_dir}/cnn.h5",
        "resnet": f"{args.models_dir}/resnet.h5",
    }

    OUT = args.outdir

    def out(name):
        return f"{OUT}/{args.prefix}{name}"

    os.makedirs(OUT, exist_ok=True)
    print(f"---- Saving to {OUT}")
    keras_model = load_model(model_files[args.model])

    input, labels = load_input(args.input, args.data_dir, args.label_column)
    input_summary = np.transpose(input.agg(func=["mean", np.sum]))
    input_summary.columns = ["mean", "sum"]

    print(f"Input of shape: {input.shape}")
    input_summary.to_csv(out("input_summary.csv"))
    print(input)
    result = predict(input, labels, keras_model)

    result["prediction"].to_csv(out("prediction.csv"), index=False)
    if result["cm"] is not None:
        result["cm"].to_csv(out("confusion_matrix.csv"))
    if result["report"] is not None:
        with open(out("report.txt"), "w") as f:
            f.write(result["report"])
    if result["metrics"] is not None:
        result["metrics"].to_csv(out("metrics.txt"))
    print(f"--- {tt.time() - start_time} seconds ---")
