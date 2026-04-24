#!/usr/bin/env python3

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Literal, NotRequired, Optional, TypedDict

# TODO: [2026-02-20 Fri] still unfinished


class FsMethod(TypedDict):
    method: NotRequired[Literal["seurat"]]
    kws: dict = field(default_factory=dict)


@dataclass
class DRSpec:
    vary: tuple[str, list[int]]
    kws: dict = field(default_factory=dict)


@dataclass
class DR:
    n_pcs: int = 50
    methods: dict[Literal["umap", "t-sne_tdr", "pacmap_tdr", "t-sne"], DRSpec] = field(
        default_factory=lambda: {"umap": ("n_neighbors", list(range(100, 5000, 500)))}
    )


@dataclass
class PermissiveFs:
    method: str = "seurat"
    kws: dict = field(
        default_factory=lambda: {
            "layer": "lshift_normalized",
            "inplace": False,
            "subset": False,
            "n_top_genes": 6000,
        }
    )


@dataclass
class DoCfs:
    cluster_exclude: list[str] | None = None  # ?this["exclude_from_clustering"]
    setup_kws: dict = field(
        default_factory=lambda: {
            "subset_hvgs": None,  # ?not test
            "leiden_kws": {},
            "pca_kws": {},
        }
    )
    cfs_kws: dict = field(default_factory=lambda: {"parallel": True})


# * QC Filtering
# ---------------------------------------------------------------------------
# mads_thresholds
# ---------------------------------------------------------------------------
# NOTE: determine these after looking at QC plots
# Map below has the following structure: COLNAME: [<OUTLIER_LOW>, <OUTLIER_HIGH>]
# A cell is determined to be an outlier if its distance from the median in terms of MADs
# falls below OUTLIER_LOW or is greater than OUTLIER_HIGH. They can be set to null not
# to do filtering in that direction


@dataclass
class MadsThresholds:
    n_genes_by_counts: list = field(default_factory=lambda: [-2, None])
    pct_counts_mito: list = field(default_factory=lambda: [None, 15])
    total_counts: list = field(default_factory=lambda: [-1, None])


@dataclass
class Normalization:
    normalize_total: dict = field(
        default_factory=lambda: {
            "target_sum": None,
            "exclude_highly_expressed": False,
        }
    )
    pooled_normalization: Optional[Any] = None
    cluster_method: str = "scran"


@dataclass
class ClusterSamples:
    features: dict[str, str | None] = field(default_factory=lambda: {"default": None})
    # The value can be a file, or list of features to use.
    # The special key `default` means to use `permissive_features` from the pipeline

    consensus_palette: str | None = None
    kws: dict = field(
        default_factory=lambda: {
            "maxK": 8,
            "reps": 1000,
            "distance": "pearson",  # 1-Pearson correlation
        }
    )


@dataclass
class DeSampleSpec:
    method: Literal["edgeR", "scVI"]
    query: str | None = None
    kws: dict = field(default_factory=dict)


class IntegrationMethod(TypedDict):
    layer: NotRequired[str]
    method: NotRequired[Literal["scVI", "harmony"]]
    kws: dict = field(default_factory=dict)


# [2026-02-20 Fri] TODO: last left to here


@dataclass
class ClusterCells:
    exclude: Optional[list[str]] = None
    method: Literal["leiden"] = "leiden"
    sweep: list[float] = field(default_factory=lambda: [0.01, 0.1, 0.5, 1])
    kws: dict = field(default_factory=lambda: {"random_state": None})


# ---------------------------------------------------------------------------
# dotplot
# ---------------------------------------------------------------------------


@dataclass
class Dotplot:
    markers: dict = field(default_factory=lambda: {})
    kws: dict = field(
        default_factory=lambda: {
            "gene_symbols": "hgnc_symbol",
            "cmap": "Greens",
            "var_group_rotation": 90,
            "layer": None,
            "log": False,
            "swap_axes": True,
        }
    )
    totals_kws: dict = field(default_factory=lambda: {"color": "green"})


# ---------------------------------------------------------------------------
# cellassign
# ---------------------------------------------------------------------------


@dataclass
class Cellassign:
    kws: dict = field(
        default_factory=lambda: {
            "layer": None,
            "size_factor_key": "lshift_size_factors",
            "train_kws": {
                "batch_size": 1024,
                "max_epochs": None,  # ?30 if test else 400
            },
        }
    )


# ---------------------------------------------------------------------------
# chosen_clusters
# ---------------------------------------------------------------------------


@dataclass
class ChosenClusters:
    leiden_res0_1: list = field(
        default_factory=lambda: ["seurat", "scVI", "leiden_res0.1"]
    )
    leiden_res0_01: list = field(
        default_factory=lambda: ["seurat", "scVI", "leiden_res0.01"]
    )


# ---------------------------------------------------------------------------
# do_de_clusters
# ---------------------------------------------------------------------------


@dataclass
class DoDeClusterScVI:
    kws: dict = field(default_factory=lambda: {"mode": "change"})


@dataclass
class DoDeClusterEdgeR:
    kws: dict = field(
        default_factory=lambda: {
            "p_value": 0.0001,
            "fc_cutoff": 1.2,
            "batch_factors": ["flowcell"],
        }
    )
    contrasts: Optional[Any] = None


@dataclass
class DoDeCluster:
    scVI: DoDeClusterScVI = field(default_factory=DoDeClusterScVI)
    edgeR: DoDeClusterEdgeR = field(default_factory=DoDeClusterEdgeR)


# ---------------------------------------------------------------------------
# enrich_clusters
# ---------------------------------------------------------------------------


@dataclass
class EnrichClusters:
    cutoff: float = 0.0001
    decouple_kws: dict = field(
        default_factory=lambda: {
            "methods": ["ulm"],
            "layer": "lshift_normalized",
        }
    )


# ---------------------------------------------------------------------------
# Top-level config
# ---------------------------------------------------------------------------


@dataclass
class Config:
    DR: DR = field(default_factory=DR)
    batch_key: str = "flowcell"
    obs_markers_annotate: list = field(
        default_factory=lambda: ["LIN28B", "GPC3", "AFP"]
    )
    mads_thresholds: MadsThresholds = field(default_factory=MadsThresholds)
    # Filtering
    normalization: Normalization = field(default_factory=Normalization)
    # Resources
    markers: Optional[Any] = None
    gene_sets: Optional[Any] = None
    gene_set_files: Any = None  # ?[f"{proj_root}/..."]
    marker_files: Any = None  # ?[f"{proj_root}/..."]

    # Feature selection
    feature_selection: dict[str, FsMethod] = field(
        default_factory=lambda: {"seurat": {"method"}}
    )
    permissive_fs: PermissiveFs = field(default_factory=PermissiveFs)

    integration: dict[str, IntegrationMethod] = field(
        default_factory=lambda: {"scVI": {"kws": {}}}
    )

    # Clustering
    do_cfs: DoCfs = field(default_factory=DoCfs)
    cluster_cells: ClusterCells = field(default_factory=ClusterCells)
    cluster_samples: ClusterSamples = field(default_factory=ClusterSamples)

    do_de_samples: dict[str, DeSampleSpec] = field(default_factory=dict)
    do_de_clusters: DoDeCluster = field(default_factory=DoDeCluster)
    enrich_clusters: EnrichClusters = field(default_factory=EnrichClusters)

    train_scvi_permissive: Any = None  # ?this["integration"]["scVI"]
    dotplot: Dotplot = field(default_factory=Dotplot)
    cellassign: Cellassign = field(default_factory=Cellassign)
    chosen_clusters: dict[str, list[str, str, str]] = field(default_factory=lambda: {})
