#!/usr/bin/env ipython

import os
from pathlib import Path

import pytest
from pyhere import here
from sarek_plus import update_gene_cache
from sarek_plus.databases import get_therapy_df
from sarek_plus.mavedb import MaveDB


def test_civic_lookup(tmp_path):
    db = tmp_path / "data.db"
    not_found = tmp_path / "not_found.txt"
    genes = {"KRAS", "EGFR", "BRAF", "RET", "dummy", "dummy2"}
    df = get_therapy_df(
        set(genes), dbs=["civic"], cache=db, not_in_db_file=not_found, lookup=True
    )
    df_prev = get_therapy_df(set(genes), cache=db, lookup=False)
    assert df.height == df_prev.height
    assert set(df["gene"]) == set(df_prev["gene"])


# def test_mavedb():


@pytest.mark.parametrize("download_path,", [(None,), ("./data/ot",)])
def test_open_targets(download_path, tmp_path):
    db = tmp_path / "data.db"
    not_found = tmp_path / "not_found.txt"
    genes = {"KRAS", "EGFR", "BRAF", "RET", "dummy", "dummy2"}
    if download_path:
        os.environ["OPEN_TARGETS_DUMP"] = download_path
    df = get_therapy_df(
        set(genes),
        dbs=["open_targets"],
        cache=db,
        not_in_db_file=not_found,
        lookup=True,
    )
    df_prev = get_therapy_df(set(genes), dbs=["open_targets"], cache=db, lookup=False)
    assert df.height == df_prev.height
    assert set(df["gene"]) == set(df_prev["gene"])


def test_database_empty(tmp_path):
    db = tmp_path / "data.db"
    not_found = tmp_path / "not_found.txt"
    genes = {"dummy", "dummy2"}
    _ = get_therapy_df(set(genes), cache=db, not_in_db_file=not_found, lookup=True)
    with open(not_found, "r") as f:
        dummies = f.read().splitlines()
        assert set(dummies) == genes


def test_update_db(tmp_path):
    db = here("tests", "sarek_plus_out", "gene_cache.db")
    os.environ["OPEN_TARGETS_DUMP"] = "./data/ot"
    import polars as pl
    from click.testing import CliRunner

    tmp_genes: Path = tmp_path / "genes1.txt"
    genes = ["KRAS", "EGFR", "BRCA", "TP53", "RET", "MET"]
    tmp_genes.write_text("\n".join(genes[:3]))
    tmp_df = pl.DataFrame({"SYMBOL": genes[3:], "AD": ["9", "2", "10"]})
    tmp_df.write_csv(tmp_path / "foobar.csv")

    runner = CliRunner()
    runner.invoke(
        update_gene_cache,
        [
            "-i",
            str(tmp_genes),
            "-i",
            str(tmp_path / "foobar.csv"),
            "-d",
            "open_targets",
            "-d",
            "civic",
            "-t",
            "-e",
            str(db),
            "-c",
            "SYMBOL",
        ],
    )
    lookup = get_therapy_df(
        set(genes), cache=db, dbs=["open_targets", "civic"], lookup=False
    )
    print(lookup)


# mdb = MaveDB("/home/shannc/Downloads/mavedb-dump.2026062418131")

# val = (
#     mdb.gene_score_set("CCR5").head(500)
#     # .select(
#     #     [
#     #         "scores.score",
#     #         "accession",
#     #         "title",
#     #         "targetGenes.HGNC",
#     #         "classification",
#     #         "hgvs.p",
#     #     ]
#     # )
# )

# from great_tables import GT, html, md

# tab = GT(val).
