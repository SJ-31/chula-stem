#!/usr/bin/env ipython

import polars as pl
import scirpy as ir
from pyhere import here

adata = ir.datasets.wu2020_3k()["airr"]
ir.pp.index_chains(adata)


vdjdb_cache = here("analyses", "data", "vdjdb.h5ad")
iedb_cache = here("analyses", "data", "iedb.h5ad")
vdjdb = ir.datasets.vdjdb(cached=True, cache_path=vdjdb_cache)

# TODO: read up on what you can use iedb for
iedb = ir.datasets.iedb(cached=True, cache_path=iedb_cache)

mcpas_path = here("analyses", "data", "McPAS-TCR.csv")
tmp = pl.read_csv(mcpas_path, infer_schema_length=None)


mcpas_cells = []
for i, row in enumerate(tmp.iter_rows(named=True)):
    cell = ir.io.AirrCell(cell_id=i)
    alpha_chain = ir.io.AirrCell.empty_chain_dict()
    beta_chain = ir.io.AirrCell.empty_chain_dict()
    alpha_chain.update(
        {
            "productive": True,
            "locus": "TRA",
            "consensus_count": 0,
            "junction": row["CDR3.alpha.nt"],
            "junction_aa": row["CDR3.alpha.aa"],
            "v_call": row["TRAV"],
            "j_call": row["TRAJ"],
            "d_call": None,
        }
    )
    beta_chain.update(
        {
            "productive": True,
            "locus": "TRB",
            "consensus_count": 0,
            "junction": row["CDR3.beta.nt"],
            "junction_aa": row["CDR3.beta.aa"],
            "v_call": row["TRBV"],
            "j_call": row["TRBJ"],
            "d_call": row["TRBD"],
        }
    )
    cell.add_chain(alpha_chain)
    cell.add_chain(beta_chain)
    mcpas_cells.append(cell)
mcpas = ir.io.from_airr_cells(mcpas_cells)
mcpas.obs = tmp.drop(
    [
        "TRBV",
        "TRBJ",
        "TRBD",
        "TRAV",
        "TRAJ",
        "CDR3.alpha.nt",
        "CDR3.alpha.aa",
        "CDR3.beta.nt",
        "CDR3.beta.aa",
    ]
).to_pandas()
ir.pp.index_chains(mcpas)
