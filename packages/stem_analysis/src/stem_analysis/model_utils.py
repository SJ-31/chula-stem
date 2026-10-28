#!/usr/bin/env ipython
#
# Utility functions intended to keep API consistent for running models
#
import argparse

import numpy as np
import pandas as pd
from sklearn.metrics import (
    f1_score,
    precision_score,
    recall_score,
)


def parse_args():
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
    parser.add_argument("--model", required=False)
    parser.add_argument(
        "--models_dir", default="models", action="store", required=False
    )
    args = parser.parse_args()
    return args


def write_input_summary(df, file, samples_x_features=True, missing_features=None):
    def count_zeros(df):
        s = df.shape[1] - df.replace({0: None}).count(axis="columns")
        s.name = "zero_null_counts"
        return s

    if samples_x_features:  # Make each row a feature
        df = np.transpose(df)
    stats = df.agg(["mean", "sum", "min", "max", "std"], axis="columns")
    zero_counts = count_zeros(df)
    stats = pd.concat([stats, zero_counts], axis="columns")
    stats["zero_null_percent"] = (stats["zero_null_counts"] / df.shape[1]) * 100
    stats["zero_null_percent"] = np.round(stats["zero_null_percent"], 2)
    stats = stats.rename_axis("feature").reset_index()
    if missing_features:
        stats["missing"] = stats["feature"].isin(missing_features)
    stats.to_csv(file, index=False)


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
