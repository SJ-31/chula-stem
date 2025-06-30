import os
from typing import Literal

import anndata as ad
from scipy import sparse

os.environ["CUDA_VISIBLE_DEVICES"] = "0"

from collections import OrderedDict

import numpy as np
import pandas as pd
import torch
from model.config import model_params as MODEL_PARAMS
from torch.utils.data import DataLoader, TensorDataset
from torch_geometric.typing import SparseTensor
from tqdm import tqdm

from utils.BulkFormer import BulkFormer

DEVICE = "cpu"

graph_path = "data/G_gtex.pt"
weights_path = "data/G_gtex_weight.pt"
gene_emb_path = "data/esm2_feature_concat.pt"


graph = torch.load(graph_path, map_location="cpu", weights_only=False)
weights = torch.load(weights_path, map_location="cpu", weights_only=False)
graph = SparseTensor(row=graph[1], col=graph[0], value=weights).t().to(DEVICE)
gene_emb = torch.load(gene_emb_path, map_location="cpu", weights_only=False)
MODEL_PARAMS["graph"] = graph
MODEL_PARAMS["gene_emb"] = gene_emb

MODEL = BulkFormer(**MODEL_PARAMS).to(DEVICE)
CKPT_MODEL = torch.load(
    "model/Bulkformer_ckpt_epoch_29.pt", weights_only=False, map_location=DEVICE
)
new_state_dict = OrderedDict()
for key, value in CKPT_MODEL.items():
    new_key = key[7:] if key.startswith("module.") else key
    new_state_dict[new_key] = value


MODEL.load_state_dict(new_state_dict)

GENE_INFO = pd.read_csv("data/bulkformer_gene_info.csv")
GENE_LIST = GENE_INFO["ensg_id"].to_list()
HIGH_VAR_GENE_IDX = torch.load("data/high_var_gene_list.pt", weights_only=False)


def extract_feature(
    expr_array,
    high_var_gene_idx,
    feature_type,
    aggregate_type,
    device,
    batch_size,
    return_expr_value=False,
    esm2_emb=None,
    valid_gene_idx=None,
):
    expr_tensor = torch.tensor(expr_array, dtype=torch.float32, device=device)
    mydataset = TensorDataset(expr_tensor)
    myloader = DataLoader(mydataset, batch_size=batch_size, shuffle=False)
    MODEL.eval()

    all_emb_list = []
    all_expr_value_list = []

    with torch.no_grad():
        if feature_type == "transcriptome_level":
            for (X,) in tqdm(myloader, total=len(myloader)):
                X = X.to(device)
                output, emb = MODEL(X, [2])
                all_expr_value_list.append(output.detach().cpu().numpy())
                emb = emb[2].detach().cpu().numpy()
                emb_valid = emb[:, high_var_gene_idx, :]

                if aggregate_type == "max":
                    final_emb = np.max(emb_valid, axis=1)
                elif aggregate_type == "mean":
                    final_emb = np.mean(emb_valid, axis=1)
                elif aggregate_type == "median":
                    final_emb = np.median(emb_valid, axis=1)
                elif aggregate_type == "all":
                    max_emb = np.max(emb_valid, axis=1)
                    mean_emb = np.mean(emb_valid, axis=1)
                    median_emb = np.median(emb_valid, axis=1)
                    final_emb = max_emb + mean_emb + median_emb

                all_emb_list.append(final_emb)
            result_emb = np.vstack(all_emb_list)
            result_emb = torch.tensor(result_emb, device="cpu", dtype=torch.float32)

        elif feature_type == "gene_level":
            for (X,) in tqdm(myloader, total=len(myloader)):
                X = X.to(device)
                output, emb = MODEL(X, [2])
                emb = emb[2].detach().cpu().numpy()
                emb_valid = emb[:, valid_gene_idx, :]
                all_emb_list.append(emb_valid)
                all_expr_value_list.append(output.detach().cpu().numpy())
            all_emb = np.vstack(all_emb_list)
            all_emb_tensor = torch.tensor(all_emb, device="cpu", dtype=torch.float32)
            esm2_emb_selected = esm2_emb[valid_gene_idx]
            esm2_emb_expanded = esm2_emb_selected.unsqueeze(0).expand(
                all_emb_tensor.shape[0], -1, -1
            )  # [B, N, D]
            esm2_emb_expanded = esm2_emb_expanded.to("cpu")

            result_emb = torch.cat([all_emb_tensor, esm2_emb_expanded], dim=-1)

    if return_expr_value:
        return np.vstack(all_expr_value_list)

    else:
        return result_emb


def main_gene_selection(X_df, gene_list):
    to_fill_columns = list(set(gene_list) - set(X_df.columns))

    padding_df = pd.DataFrame(
        np.full((X_df.shape[0], len(to_fill_columns)), -10),
        columns=to_fill_columns,
        index=X_df.index,
    )

    X_df = pd.DataFrame(
        np.concatenate([df.values for df in [X_df, padding_df]], axis=1),
        index=X_df.index,
        columns=list(X_df.columns) + list(padding_df.columns),
    )
    X_df = X_df[gene_list]

    var = pd.DataFrame(index=X_df.columns)
    var["mask"] = [1 if i in to_fill_columns else 0 for i in list(var.index)]
    return X_df, to_fill_columns, var


def demo():
    demo_df = pd.read_csv("data/demo.csv")
    n_samples = 3

    input_df, to_fill_columns, var = main_gene_selection(
        X_df=demo_df, gene_list=GENE_LIST
    )

    var.reset_index(inplace=True)
    valid_gene_idx = list(var[var["mask"] == 0].index)

    # Extract transcritome-level embedding
    result = extract_feature(
        expr_array=input_df.values[:n_samples],
        high_var_gene_idx=HIGH_VAR_GENE_IDX,
        feature_type="transcriptome_level",
        aggregate_type="max",
        device=DEVICE,
        batch_size=4,
        return_expr_value=False,
        esm2_emb=MODEL_PARAMS["gene_emb"],
        valid_gene_idx=valid_gene_idx,
    )
    print(result.shape)

    # Extract gene-level embedding
    result = extract_feature(
        expr_array=input_df.values[:n_samples],
        high_var_gene_idx=HIGH_VAR_GENE_IDX,
        feature_type="gene_level",
        aggregate_type="all",
        device=DEVICE,
        batch_size=4,
        return_expr_value=False,
        esm2_emb=MODEL_PARAMS["gene_emb"],
        valid_gene_idx=valid_gene_idx,
    )

    # Extract expression values
    result = extract_feature(
        expr_array=input_df.values[:n_samples],
        high_var_gene_idx=HIGH_VAR_GENE_IDX,
        feature_type="transcriptome_level",
        aggregate_type="all",
        device=DEVICE,
        batch_size=4,
        return_expr_value=True,
        esm2_emb=MODEL_PARAMS["gene_emb"],
        valid_gene_idx=valid_gene_idx,
    )
    print(result.shape)


def log_tpm(adata: ad.AnnData) -> np.ndarray:
    expr = adata.X if not sparse.issparse(adata.X) else adata.X.toarray() + 1
    lengths = adata.var["SEQLENGTH"]
    numer = np.log(expr) - np.reshape(np.log(lengths), (1, -1))
    denom = np.log(np.nansum(np.exp(numer), axis=1)).reshape(-1, 1)
    tpm = np.exp(numer - denom + np.log(1e6))
    tpm = np.nan_to_num(tpm, neginf=0)
    return tpm


def extract_helper(
    adata: ad.AnnData,
    feature_type: Literal["transcriptome_level", "gene_level"] = "transcriptome_level",
    aggregate_type: Literal["max", "mean", "median", "all"] = "max",
):
    expr = log_tpm(adata)
    df = pd.DataFrame(expr, columns=adata.var["GENEID"])
    input_df, _, var = main_gene_selection(X_df=df, gene_list=GENE_LIST)
    var.reset_index(inplace=True)
    valid_gene_idx = list(var[var["mask"] == 0].index)
    result: np.ndarray = extract_feature(
        expr_array=input_df.values,
        high_var_gene_idx=HIGH_VAR_GENE_IDX,
        feature_type=feature_type,
        aggregate_type=aggregate_type,
        device=DEVICE,
        batch_size=4,
        return_expr_value=False,
        esm2_emb=MODEL_PARAMS["gene_emb"],
        valid_gene_idx=valid_gene_idx,
    )
    if isinstance(result, torch.Tensor):
        result = result.detach().numpy()
    if feature_type != "transcriptome_level":
        var = pd.DataFrame({"GENEID": GENE_LIST})
        filtered = adata.var.loc[adata.var["GENEID"].isin(GENE_LIST), :]
        var = var.merge(filtered, how="left", on="GENEID")
    else:
        var = None
    new = ad.AnnData(X=sparse.csc_array(result), obs=adata.obs, var=var)
    return new


def parse_args():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input")
    parser.add_argument("-o", "--output")
    parser.add_argument("-a", "--aggregate_type", default="max", action="store")
    parser.add_argument(
        "-f",
        "--feature_type",
        default="transcriptome_level",
        help="What feature type to extract",
        action="store",
    )
    parser.add_argument("-d", "--demo", action="store_true")
    args = vars(parser.parse_args())
    return args


if __name__ == "__main__":
    args = parse_args()
    if args["demo"]:
        demo()
    else:
        adata = ad.read_h5ad(args["input"])
        embedded: ad.AnnData = extract_helper(
            adata,
            feature_type=args["feature_type"],
            aggregate_type=args["aggregate_type"],
        )
        embedded.write_h5ad(args["output"])
