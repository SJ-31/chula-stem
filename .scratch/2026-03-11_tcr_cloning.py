#!/usr/bin/env ipython

from pathlib import Path

import mudata as mu
import pandas as pd
import polars as pl
import polars.selectors as cs
from great_tables import GT
from pyhere import here

outdir: Path = here("analyses", "output", "pdac_tcr", "2026-01-22")

airr = mu.read_h5mu(outdir / "combined.h5mu")["airr"]

top_clones = (
    pl.read_csv(outdir / "top_clones.csv")
    .unique(["Sample_Name", "clone_id"])
    .filter(pl.col("rank") <= 3)
    .unique(["Sample_Name", "rank"])
    .sort(["Sample_Name", "count"], descending=True)
)

count_tab = GT(
    airr.obs.groupby("Sample_Name").agg({"clone_id": "nunique"}).reset_index()
).cols_label({"clone_id": "Number of clones"})
count_tab.save(outdir / "count_table.png", scale=5)

top_tab = (
    GT(top_clones)
    .cols_hide(columns=cs.contains("_2"))
    .tab_stub(groupname_col="Sample_Name")
    .tab_stubhead(label="Sample")
    .tab_options(row_group_as_column=True)
)
top_tab.show()
top_tab.save(outdir / "top_table.png", scale=5)
