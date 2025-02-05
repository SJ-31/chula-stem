#!/usr/bin/env python
import argparse
import os

import numpy as np
import pandas as pd
import torch
from sklearn.metrics import (
    classification_report,
    confusion_matrix,
    f1_score,
    precision_score,
    recall_score,
)

import model_utils as mu
from dataset import MyDataSet
from model import Pathway_Guided_Transformer
from utils import TCGA_labelDict, evaluate

LABEL_DICT = TCGA_labelDict
LABEL_DICT_REV = {v: k for k, v in LABEL_DICT.items()}

DEVICE = torch.device("cuda:2" if torch.cuda.is_available() else "cpu")


def reshape_data(data):
    return data.reshape(data.shape[0], 1, data.shape[1])


def raw_to_kegg(df, data_dir: str):
    kegg_info = pd.read_csv(f"{data_dir}/KEGG_Pathway_information.csv", header=0)
    all_gene_list = []
    for gene_set in kegg_info["gene_name"]:
        all_gene_list.extend(gene_set.split("/"))

    gene_num = 1
    kegg_data_dict = dict()
    for gene in all_gene_list:
        kegg_data_dict["gene" + str(gene_num)] = df[gene]
        gene_num = gene_num + 1
    converted = pd.DataFrame(kegg_data_dict)
    return converted


def load_input(file: str, label_column: str, data_dir: str):
    """`file` is either a gene x sample matrix, or a sample x gene matrix
        the latter case is if labels are known (in a column called "tumor_type")
        and the file is provided for testing purposes.

    entries are gene expression values in log2 fpkm

    :returns: a sample x gene matrix (which the model expects),
    """
    df = pd.read_csv(file, index_col=0)
    if label_column in df.columns:
        labels: pd.Series = df[label_column]
        df = np.transpose(df.drop([label_column], axis=1))
    else:
        labels = pd.Series(None, index=df.columns, name=label_column)

    rows = []
    with open(f"{data_dir}/features.txt", "r") as f:
        target_genes = f.read().splitlines()

    missing = set()
    for t in target_genes:
        try:
            rows.append(df.loc[[t],])
        except KeyError:
            missing.add(t)
            print(f"WARNNG: feature {t} missing")
            rows.append(pd.DataFrame(None, columns=df.columns, index=[t]))
    selected = pd.concat(rows, axis=0)
    transposed: pd.DataFrame = np.transpose(selected)
    converted = raw_to_kegg(transposed, data_dir)
    return converted, labels, missing, transposed


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
    input: pd.DataFrame,
    labels: pd.Series,
    path_weight,
    batch_size=64,
):
    samples = input.index
    reshaped = reshape_data(input.values)  # Shape of (n_samples, 1, n_genes)
    # TODO: <2025-01-31 Fri> Need to make a version of this that doesn't have labels
    encoded_labels = [TCGA_labelDict[i] for i in labels]
    to_test = MyDataSet(reshaped, encoded_labels)
    test_loader = torch.utils.data.DataLoader(to_test, batch_size=batch_size)
    model.load_state_dict(torch.load(path_weight, map_location=DEVICE))
    model.eval()
    criterion = torch.nn.CrossEntropyLoss()
    test_loss, test_acc, pre_label = evaluate(
        model=model,
        data_loader=test_loader,
        num=len(to_test),
        criterion=criterion,
        device=DEVICE,
    )
    predictions = [LABEL_DICT_REV[i] for i in pre_label]
    with_pred = pd.DataFrame({"prediction": predictions, "sample": samples})
    result = {"prediction": with_pred, "cm": None, "report": None, "metrics": None}
    if labels.any():
        metrics = sk_report_performance(labels, predictions)
        print(metrics)
        result["metrics"] = metrics
        result["report"] = classification_report(labels, predictions)
        sorted_labels = list(set(predictions) | set(labels))
        sorted_labels.sort()
        result["cm"] = pd.DataFrame(
            confusion_matrix(labels, predictions, labels=sorted_labels),
            index=sorted_labels,
            columns=sorted_labels,
        )
    return result


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input",
        default="",
        required=True,
        help="csv file containing expression data",
        action="store",
    )
    parser.add_argument(
        "-b", "--batch_size", default=64, help="Model batch size", action="store"
    )
    parser.add_argument(
        "-w",
        "--weights",
        required=False,
        default="best_weight.pth",
        help="Path to weights",
    )
    parser.add_argument("--model", default="", help="", action="store")
    parser.add_argument("--data_dir", default="data", action="store")
    parser.add_argument("-p", "--prefix", default="cupai", help="Prefix for output")
    parser.add_argument("-o", "--outdir", default=".", help="Output directory")
    parser.add_argument(
        "-l",
        "--label_column",
        default="label",
        help="For test datasets (sample x genes), the column indicating the label of the sample",
        action="store",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    pathway_df = pd.read_csv("RNAseq/KEGG_Pathway_information.csv", header=0)
    pathway_num = list(pathway_df["count"])
    model = Pathway_Guided_Transformer(
        num_classes=32,
        pathway_number=pathway_num,
        dim=512,
        depth=6,
        heads=8,
        mlp_dim=1024,
        dropout=0.1,
        emb_dropout=0.1,
    ).to(DEVICE)
    model = Pathway_Guided_Transformer(
        num_classes=32,
        pathway_number=pathway_num,
        dim=512,
        depth=6,
        heads=8,
        mlp_dim=1024,
        dropout=0.1,
        emb_dropout=0.1,
    ).to(DEVICE)

    def out(name):
        return f"{args.outdir}/{args.prefix}{name}"

    os.makedirs(args.outdir, exist_ok=True)
    data, labels, missing, before = load_input(
        args.input, args.label_column, args.data_dir
    )
    mu.write_input_summary(data, out("input_summary_kegg.csv"))
    mu.write_input_summary(
        before, out("input_summary_genes.csv"), missing_features=missing
    )

    result = predict(data, labels, args.weights, args.batch_size)

    result["prediction"].to_csv(out("prediction.csv"), index=False)
    if result["cm"] is not None:
        result["cm"].to_csv(out("confusion_matrix.csv"))
    if result["metrics"] is not None:
        result["metrics"].to_csv(out("metrics.csv"))
    if result["report"] is not None:
        with open(out("report.txt"), "w") as f:
            f.write(result["report"])
