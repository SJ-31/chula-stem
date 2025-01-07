#!/usr/bin/env ipython

import pytest
from chula_stem.databases import get_therapy_df


def test_database_lookup():
    civic = "./civic.json"
    pd = "./pandrugs2.json"
    not_found = "./data/not_in_db.txt"
    with open("./data/sample_genes.txt", "r") as f:
        gene_list = f.read().splitlines()
    get_therapy_df(
        gene_list, civic_cache=civic, pandrugs2_cache=pd, not_in_db_file=not_found
    )
