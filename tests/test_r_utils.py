#!/usr/bin/env ipython

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scanpy as sc
from chula_stem.r_utils import edgeR_wrapper


@pytest.fixture
def toy_data() -> ad.AnnData:
    adata: ad.AnnData = sc.datasets.pbmc68k_reduced()
    summed = sc.get.aggregate(adata, by=["bulk_labels", "phase"], func="sum")
    summed.X = np.round(np.abs(summed.layers["sum"]), 0)
    return summed


def test_ovr(toy_data):
    result = edgeR_wrapper(toy_data, group="phase", id_col=None)


def test_contrasts(toy_data):
    result = edgeR_wrapper(
        toy_data,
        group="phase",
        id_col=None,
        contrasts={"g1_vs_s": "phaseG1-phaseS", "g2m_vs_g1": "phaseG2M-phaseG1"},
    )
    print(result)


def test_avr(toy_data):
    result = edgeR_wrapper(
        toy_data,
        group="phase",
        id_col=None,
        ovr=False,
        avr={
            "g1_vs_all": {"phase": ["G1"]},
            "s_g1_vs_all": {"phase": ["G1", "S"]},
        },
    )
    print(result)
