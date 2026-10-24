#!/usr/bin/env ipython

from pathlib import Path

import anndata as ad
import chula_stem.pdac_subtyping as pds
import pandas as pd
from chula_stem.sc_rnaseq import annotate_adata_vars
from pyhere import here
from scipy import sparse

ad.settings.allow_write_nullable_strings = True

tx2gene = pd.read_csv(here("analyses", "data", "tx2gene.tsv"), sep="\t")
workdir = here("analyses", "pdac_subtyping")
pdata = here("analyses", "data_all", "public_data")

manifest = pd.read_csv(here(workdir, "manifest.csv"))

adata_path: Path = here(workdir, "cohort.h5ad")

if not adata_path.exists():
    adata = pds.read_manifest(manifest, tx2gene, purecn_prefix="")
    adata = annotate_adata_vars(
        adata,
        savepath=here("analyses", "data", "ensembl_gene_data.csv"),
        new_index_col="hgnc_symbol",
        new_index_name="gene",
    )
    adata.X = sparse.csc_array(adata.X)
    adata.obs.loc[:, "sample"] = adata.obs_names
    adata.write_h5ad(adata_path)
else:
    adata = ad.read_h5ad(adata_path)

pds.add_kras_imbalance(manifest, adata, purecn_prefix="")
pds.add_waddell_classification(manifest, adata)
model = pds.MoffittBasal(threshold=0.44)  # From running last time
model.fit(adata=adata)
model.predict(inplace=True)
print(adata.obs)
adata.write_h5ad(adata_path)
for key in ("KRAS_imbalance_subtype", "moffitt_is_basal", "waddell_subtype"):
    adata.obs.loc[:, ["sample", key]].rename({key: "label"}, axis="columns").to_csv(
        workdir / f"{key}.tsv", index=False, sep="\t", na_rep="NA"
    )
pds.save_purecn_summary(adata.uns["purecn_summary"], workdir / "purecn_summary.png")
