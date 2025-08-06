#!/usr/bin/env python
import re
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_dir", default=".", action="store")
    parser.add_argument(
        "-s", "--suffix", default="Msisensor_summary.tsv", action="store"
    )
    parser.add_argument("-o", "--output")
    args = vars(parser.parse_args())  # convert to dict
    return args


if __name__ == "__main__":
    args = parse_args()
    msi_files = Path(args["input_dir"]).rglob(f"*{args["suffix"]}")
    dfs = []
    for file in msi_files:
        sample = re.findall(f"[0-9]+-(.*)-{args['suffix']}", file.name)[0]
        df = pd.read_csv(file, sep="\t").assign(sample=sample)
        if df.shape[0] == 0:
            df = pd.DataFrame({k: np.nan for k in df.columns}, index=[0])
        dfs.append(df)
    pd.concat(dfs).to_csv(args["output"], index=False, sep="\t", na_rep="NA")
