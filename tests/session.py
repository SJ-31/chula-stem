#!/usr/bin/env ipython

from importlib import util
from pathlib import Path
from typing import Callable, Literal

import anndata as ad
import mudata as md
import numpy as np
import polars as pl
import polars.selectors as cs
import pydna.tm as tm
import scirpy as ir
from great_tables import GT
from loguru import logger
from pydna.design import Amplicon, pcr, primer_design
from pydna.dseqrecord import Dseqrecord
from pydna.primer import Primer
from pyhere import here
from skbio import DNA
from skbio.alignment import PairAlignPath, pair_align


def load_file_as_module(name, location):
    spec = util.spec_from_file_location(name, location)
    module = util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


mod = load_file_as_module("primers", here("snakemake", "tcr", "scripts", "primers.py"))

# %%

combined = md.read_h5mu(
    here("analyses", "output", "pdac_tcr", "2026-01-22", "combined.h5mu")
)


# NOTE: bd rhapsody data doesn't seem to give the indices for cdr3/junction, so
# we will get this manually by alignment
def locate_cdr3(
    cdr3: str,
    full: str,
    j_end: int | None = None,
    v_start: int | None = None,
    exact_match: bool = True,
) -> tuple[int, int]:
    """Helper function to locate CDR3 within full-length AIR sequence
    TODO: you need this to actually use the biology. It will fail if there
    is a subsequence of v gene that happens to be identical to cdr3


    Parameters
    ----------
    exact_match : bool

    Notes
    -----
    We expect that the cdr3

    """
    try:
        cdr3_start: int = full.index(cdr3)
        cdr3_end = cdr3_start + len(cdr3)
    except ValueError:
        if exact_match:
            raise ValueError("Exact match for cdr3 sequence in `full` not found")
        logger.warning("Exact match for cdr3 not found. Resorting to alignment")
        res = pair_align(cdr3, full, mode="global", free_ends=(True, False))
        path: PairAlignPath = res.paths[0]
        # "global" is guaranteed to have at least one path
        coords: np.ndarray = path.to_coordinates()
        cdr3_start, cdr3_end = coords[1].tolist()[1:3]
    # if v_start and cdr3_start <= v_start:
    #     logger.info(v_start, cdr3_start)
    #     raise ValueError("cdr3 sequence start invalid: before start of V gene")
    # if j_end and cdr3_end >= j_end:
    #     raise ValueError("cdr3 sequence end invalid: after end of J gene")
    return cdr3_start, cdr3_end


wanted_clones = [["136", "PDAC83"], ["354", "PDAC82"]]
df = mod.filter_wanted_clones(
    combined["airr"],
    {"wanted_clones": wanted_clones},
    extras=["fwr1", "fwr2", "fwr3", "fwr4", "cdr2", "cdr1"],
)

# * Testzone
# [2026-01-29 Thu] TODO: All you want is something that can index the cdr3...
ww = next(df.iter_rows(named=True))
cdr3 = ww["VJ_1_cdr3"]
lcdr3 = len(cdr3)
v_start = ww["VJ_1_v_sequence_start"]
v_end = ww["VJ_1_v_sequence_end"]
cdr3_start = ww["VJ_1_cdr3_start"]
cdr3_end = ww["VJ_1_cdr3_end"]
v_end = ww["VJ_1_v_sequence_end"]
j_end = ww["VJ_1_j_sequence_end"]
full = ww["VJ_1_sequence"]
v_seq = full[v_start:v_end]
up_to = cdr3_end - v_start


forward = full[v_start : v_start + 30]
reverse = str(DNA(full[cdr3_start:cdr3_end]).reverse_complement())


cur, dct = mod.format_one_sample(ww, cfg={"end_gene": "j"})


res = mod.assemble_one_construct(
    ww,
    ["VJ_1", "VDJ_1"],
    {
        "include_leader": True,
        "include_c": True,
        "sequences": {
            "linker": "TTA",
            "flanking": {"five_prime": "TACGCCGCAT", "three_prime": "CGCGAGCTA"},
            "c_gene": {
                "TRB": "/home/shannc/Bio_SDD/chula-stem/analyses/pdac_tcr/construct_sequences/murinized_trbc.fasta",
                "TRA": "/home/shannc/Bio_SDD/chula-stem/analyses/pdac_tcr/construct_sequences/murinized_trac.fasta",
            },
        },
    },
)
names = res["named"]
validation = res["validation"]
ax, _ = res["plot"].plot(figure_width=5)
ax.figure.show()
cat = DNA.concat(names.values())

# BUG: resolve C gene stop codon issue
