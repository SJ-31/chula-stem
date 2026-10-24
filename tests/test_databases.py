#!/usr/bin/env ipython

from chula_stem.databases import get_therapy_df


def test_database_lookup(tmp_path):
    db = tmp_path / "data.db"
    not_found = tmp_path / "not_found.txt"
    with open("./data/sample_genes.txt", "r") as f:
        gene_list = f.read().splitlines()
    df = get_therapy_df(set(gene_list), cache=db, not_in_db_file=not_found, lookup=True)
    print(df)
    print(df.schema)
