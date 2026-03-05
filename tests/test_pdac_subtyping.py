#!/usr/bin/env ipython

from pathlib import Path

import anndata as ad
import chula_stem.pdac_subtyping as pds
import h5py
import pandas as pd
from chula_stem.r_utils import tximport
from pyhere import here

tx2gene = pd.read_csv(here("analyses", "data", "tx2gene.tsv"), sep="\t")
workdir = here("analyses", "pdac_subtyping")

manifest = pd.read_csv(here(workdir, "manifest.csv"))

adata_path: Path = here(workdir, "pdac_cohort.h5ad")

if not adata_path.exists():
    adata = pds.read_manifest(manifest, tx2gene)
    adata.write_h5ad(adata_path)
else:
    adata = ad.read_h5ad(adata_path)
