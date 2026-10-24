#!/usr/bin/env ipython

from chula_stem.databases import get_therapy_df


def test_database_lookup(tmp_path):
    db = tmp_path / "data.db"
    not_found = tmp_path / "not_found.txt"
    genes = {"KRAS", "EGFR", "BRAF", "RET", "dummy", "dummy2"}
    df = get_therapy_df(set(genes), cache=db, not_in_db_file=not_found, lookup=True)
    df_prev = get_therapy_df(set(genes), cache=db, lookup=False)
    assert df.height == df_prev.height
    assert set(df["gene"]) == set(df_prev["gene"])
    print(df)
    print(df.schema)


def test_database_empty(tmp_path):
    db = tmp_path / "data.db"
    not_found = tmp_path / "not_found.txt"
    genes = {"dummy", "dummy2"}
    _ = get_therapy_df(set(genes), cache=db, not_in_db_file=not_found, lookup=True)
    with open(not_found, "r") as f:
        dummies = f.read().splitlines()
        assert set(dummies) == genes
