#!/usr/bin/env python3

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


def xarray_if_sparse(
    x: ad.AnnData, copy: bool = True, dtype: np.dtype = np.float32
) -> np.ndarray:
    was_sparse: bool = sparse.issparse(x.X)
    if was_sparse:
        arr = x.X.toarray()
    elif copy:
        arr = x.X.copy()
    else:
        arr = x.X
    return arr.astype(dtype)


def enrichment_set_cover(
    pathways: dict[str, set],
    p_values: dict[str, float],
    coverage: float = 0.8,
) -> tuple[set, float]:
    """
    Apply set cover to remove redundant enriched pathways

    Parameters
    ----------
    pathways : dict[str, set]
        Dictionary mapping names of enriched pathways to their member genes
    p_values : dict[str, float]
        Dict mapping pathways to their (adjusted) p values

    Notes
    -----
    The algorithm applies the set cover algorithm to prioritize pathways
    with the lowest p-values

    Stoney, R.A., Schwartz, JM., Robertson, D.L. et al. Using set theory to reduce redundancy in pathway sets. BMC Bioinformatics 19, 386 (2018). https://doi.org/10.1186/s12859-018-2355-3
    """
    universe: set = {s for st in pathways.values() for s in st}
    uncovered: set = universe.copy()
    covered: set = set()
    result: set = set()

    def score_set(name: str, s: set) -> int:
        if len(s & uncovered) != 0:
            return 1 - p_values[name]
        return 0

    while (len(covered) / len(universe)) < coverage:
        scored = [(score_set(name, s), name) for name, s in pathways.items()]
        top_score, tname = sorted(scored, reverse=True, key=lambda x: x[0])[1]
        in_top = pathways[tname]
        covered |= in_top
        uncovered -= in_top
        del pathways[tname]
        result.add(tname)
    return result, len(covered) / len(universe)


def set_cover_on_dc_results(
    df: pd.DataFrame, net: pd.DataFrame, coverage: float
) -> pd.DataFrame:
    """Apply `enrichment_set_cover` to decoupler results after ranking
    with `dc.tl.rankby_group`
    """
    net_mapping: dict[str, set] = (
        net.groupby("source").agg({"target": set}).to_dict()["target"]
    )
    result = []
    for name, group in df.groupby(["group", "reference"]):
        pwy = {n: s for n, s in net_mapping.items() if n in group["name"].values}
        pvals: dict = (
            group.loc[:, ["name", "padj"]]
            .drop_duplicates("name")
            .set_index("name")
            .to_dict()["padj"]
        )
        covering, score = enrichment_set_cover(
            pathways=pwy, p_values=pvals, coverage=coverage
        )
        result.append(group.loc[group["name"].isin(covering), :])
    return pd.concat(result)
