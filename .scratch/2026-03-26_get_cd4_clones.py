#!/usr/bin/env python3

import muon as mu
from pyhere import here

data = mu.read_h5mu(
    here("analyses", "output", "pdac_tcr", "2026-01-22", "combined.h5mu")
)

airr = data.mod["airr"]

# Get two of the top clones CD4 memory from PDAC83
rna_obs = data.mod["rna"].obs

cell_type_map = rna_obs.loc[:, ["Cell_Type_Experimental"]]

airr_obs = airr.obs.merge(cell_type_map, left_index=True, right_index=True)

e3 = airr_obs.loc[airr.obs["Sample_Name"] == "PDAC83", :]

print(e3.Cell_Type_Experimental.value_counts())

print(
    e3.loc[:, ["clone_id", "clone_id_size", "Cell_Type_Experimental"]]
    .sort_values("clone_id_size")
    .drop_duplicates("clone_id")
    .tail()
)
