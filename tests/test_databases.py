#!/usr/bin/env ipython

import pytest
from chula_stem.databases import get_therapy_df


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


@pytest.mark.parametrize("download_path,", [None, "./data/ot"])
def test_open_targets(download_path, tmp_path):
    db = tmp_path / "data.db"
    not_found = tmp_path / "not_found.txt"
    genes = {"KRAS", "EGFR", "BRAF", "RET", "dummy", "dummy2"}
    df = get_therapy_df(
        set(genes),
        dbs=["open_targets"],
        cache=db,
        not_in_db_file=not_found,
        lookup=True,
        db_kws={"open_targets": {"data_dir": download_path}},
    )
    print(df)
    df_prev = get_therapy_df(set(genes), cache=db, lookup=False)
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
