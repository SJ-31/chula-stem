#!/usr/bin/env ipython

import polars as pl
import polars.selectors as cs
from pyhere import here

wanted_cells = {
    # CD8
    "Exhausted CD8+ T cell",
    "Effector CD8+ T cell",
    "Cytotoxic CD8+ T cell",
    "Memory CD8+ T cell",
    "Activated CD8+ T cell",
    # CD4
    "Memory CD4+ T cell",
    "Naive CD4+ T cell",
    "Activated CD4+ T cell",
    "CD4+ T cell",
    # B cell
    "Activated B cell",
    "Memory B cell",
    "Naive B cell",
    # Misc
    "Natural killer T (NKT) cell",
    "Gamma delta T cell",
    "Non-classical monocyte",
    "Classical monocyte",
}
to_replace = {
    "Naive CD4 T cell": "Naive CD4+ T cell",
    "Gamma delta(Î³Î´) T cell": "Gamma delta T cell",
    "Effector CD8+ memory T (Tem) cell": "Memory CD8+ T cell",
}


# * Format df

cellmarker = (
    pl.read_csv(here("analyses", "data", "CellMarker2_human.csv"))
    .filter((pl.col("cell_type") == "Normal cell") & pl.col("Symbol").is_not_null())
    .with_columns(pl.col("cell_name").replace(to_replace))
)

filtered = cellmarker.filter(
    pl.col("cell_name").is_in(wanted_cells) & pl.col("Symbol").is_not_null()
)


binary = (
    filtered.select(["cell_name", "Symbol"])
    .unique()
    .with_columns(pl.lit(1).alias("value"))
    .pivot(on="cell_name", index="Symbol", values="value")
    .with_columns(cs.numeric().fill_null(0))
)

meta_output = here("analyses", "pdac_tcr", "cellassign_markers_meta.csv")
binary_file = here("analyses", "pdac_tcr", "cellassign_markers.csv")

filtered.write_csv(meta_output)
binary.write_csv(binary_file)
