#!/usr/bin/env ipython

from collections.abc import Sequence

import polars as pl
import scirpy as ir
import yaml
from pyhere import here

with open(here("analyses", "pdac_tcr", "env.yaml"), "r") as f:
    config = yaml.safe_load(f)

source("/home/shannc/Bio_SDD/chula-stem/snakemake/tcr/prepare_data.py")
source("/home/shannc/Bio_SDD/chula-stem/snakemake/tcr/alignment.py")
source("/home/shannc/Bio_SDD/chula-stem/snakemake/tcr/plotting.py")
source("/home/shannc/Bio_SDD/chula-stem/snakemake/tcr/analyses.py")


def find_covered_ranges(seq: Sequence):
    ranges = []
    current = []
    for i, char in enumerate(seq):
        if char is not None:
            current.append(i)
        if (char is None or i == len(seq) - 1) and current:
            current.append(i - 1)
            ranges.append((min(current), max(current)))
            current = []
    return ranges


mdata = md.read_h5mu(
    "/home/shannc/Bio_SDD/chula-stem/analyses/output/pdac_tcr/no_date/combined.h5mu"
)

airr = mdata["airr"]
sample = "PDAC82"

ori = airr[airr.obs["Sample_Name"] == sample, :]


_, filtered = maybe_filter_by_rank(ori, config["mixtcrpred"])


cols = ["sequence", "cdr3", "cdr1", "cdr2", "fwr1", "fwr2", "fwr3", "fwr4"] + [
    s
    for seq in [
        [f"{k}_sequence_start", f"{k}_sequence_end", f"{k}_call"]
        for k in ("v", "d", "j")
    ]
    for s in seq
]
seqs = pl.from_pandas(ir.get.airr(filtered, ["sequence_id"] + cols + ["c_call"]))


nucleotide_colors = {"A": "red", "T": "blue", "C": "green", "G": "yellow"}

# [2025-10-15 Wed] You can do this but tbh it's not worth the effort
row = seqs.row(0, named=True)
patterns = {k: row.get(f"VJ_1_{k}") for k in SEQ_COLS}
test = global_local_alignment(patterns, patterns.pop("sequence"))
test["sequence"] = test[""]
del test[""]

align_df = pl.DataFrame({k: list(v) for k, v in test.items()}).with_columns(
    pl.all().replace("-", None)
)
font_size = 12
wrap_length = 400
seq = align_df["sequence"]
seqlen = len(seq)
align_df = align_df.drop("sequence")
nrows = seqlen // wrap_length
fig, ax = plt.subplots(nrows=max(nrows + 1, 2))
ax[0].set_xlim(0, wrap_length)
# ax[0].axis("off")
ax[0].set_ylim(0, 5)
tracked = []
row_tracker, counter = 0, 0
separation = 3
for i, (char, row) in enumerate(zip(seq, align_df.iter_rows(named=True))):
    incr = 1
    if counter >= wrap_length:
        row_tracker += 1
        ax[row_tracker].set_xlim(0, wrap_length)
        # ax[row_tracker].axis("off")
        ax[row_tracker].set_ylim(0, 5)
        counter = 0
    ax[row_tracker].text(
        counter,
        0,
        char,
        fontsize=font_size,
        backgroundcolor=nucleotide_colors.get(char),
    )
    for k, v in row.items():
        if v is not None:
            ax[row_tracker].text(counter, incr, v, fontsize=font_size)
            incr += 1
    tracked.append(incr)
    counter += 1

for col in align_df.columns:
    ranges = find_covered_ranges(align_df[col])
    for i, rnge in enumerate(ranges):
        ax[row_tracker].axhline(i, rnge[0], rnge[1], label=col, lw=0.8)

fig.subplots_adjust(hspace=0.6)
fig.set_figwidth(100)
fig.set_dpi(500)
fig.savefig("/home/shannc/Downloads/dummy.png")
