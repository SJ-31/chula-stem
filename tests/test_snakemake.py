#!/usr/bin/env ipython


import mudata as md
import numpy as np
import polars as pl
import pytest
import scirpy as ir
from pyhere import here


def load_file_as_module(name, location):
    from importlib import util

    spec = util.spec_from_file_location(name, location)
    module = util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# * TCR pipeline


tcr_primers = load_file_as_module(
    "tcr", here("snakemake", "tcr", "scripts", "primers.py")
)

tcr_cons = load_file_as_module(
    "tcr_cons", here("snakemake", "tcr", "scripts", "construct.py")
)


@pytest.fixture
def airr_data_test():
    combined = md.read_h5mu(here("tests", "data", "airr_test.h5mu"))
    return combined["airr"]


def test_get_c(airr_data_test):
    chain = "VJ_1"
    df = pl.from_pandas(
        ir.get.airr(airr_data_test, tcr_primers.AIRR_WANTED_COLS, chain=chain),
        include_index=True,
    )
    dct = next(df.iter_rows(named=True))
    full = dct[f"{chain}_sequence"]
    j_end = f"{chain}_j_sequence_end"
    j_start = f"{chain}_j_sequence_start"
    old_j = full[dct[j_start] : dct[j_end]]
    old_c = full[dct[j_end] :]
    cdata = tcr_primers.ChainData(chain, dct)
    new_c = tcr_primers.infer_air_endpoint(cdata, c_try_match_allele=True)
    print(f"Old: {old_c}\nNew: {new_c}")
    print(f"Old full: {full}")
    new_full = tcr_primers.add_new_c(cdata, str(new_c))
    print(f"New full: {new_full}")
    assert old_j == new_full[dct[j_start] : dct[j_end]]


def create_dummy_construct(
    five_prime, three_prime, leader, c_gene, linker
) -> tuple[dict, dict, str]:
    """
    Helper function to create a random TCR construct from the configuration

    Returns the tcr construct and the airr data dictionary with the fwr, cdr regions
    populated. Genes are filled in randomly (do not bother testing the visualization)
    """
    cfg = {"sequences": {"flanking": {}, "c_gene": {}, "leader": {}}}
    result = {
        "VJ_1_v_call": "TRA*434",
        "VDJ_1_v_call": "TRB*423",
        "VDJ_1_j_call": "TRBJ1*01",
        "VJ_1_j_call": "TRAJ1*01",
    }
    seq = []
    if five_prime:
        seq.append(five_prime)
        cfg["sequences"]["flanking"]["five_prime"] = five_prime
    if c_gene:
        cfg["sequences"]["c_gene"]["TRA"] = c_gene
        cfg["sequences"]["c_gene"]["TRB"] = c_gene
        cfg["include_c"] = True
    else:
        cfg["include_c"] = False
    if linker:
        cfg["sequences"]["linker"] = linker
    if leader:
        cfg["sequences"]["leader"]["TRA"] = leader
        cfg["sequences"]["leader"]["TRB"] = leader
        cfg["include_leader"] = True
    else:
        cfg["include_leader"] = False
    regions = ("fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3", "fwr4")
    for i, chain in enumerate(("VJ_1", "VDJ_1")):
        if leader:
            seq.append(leader)
        for region in regions:
            cur = "".join(
                list(np.random.choice(["A", "T", "C", "G"], size=5, replace=True))
            )
            result[f"{chain}_{region}"] = cur
            seq.append(cur)
        if c_gene and (not c_gene.endswith("TAA") or i == 1):
            seq.append(c_gene)
        elif c_gene:
            seq.append(c_gene[:-3])
        if linker and i != 1:
            seq.append(linker)
    if three_prime:
        cfg["sequences"]["flanking"]["three_prime"] = three_prime
        seq.append(three_prime)
    return cfg, result, "".join(seq)


@pytest.mark.parametrize(
    "five_prime,three_prime,leader,c_gene,linker",
    [
        (None, None, None, None, None),
        ("ATCGG", None, "GCTAT", "TAGCA", None),
        (None, "CGATC", "TTGAC", None, "ACGTA"),
        ("GGCAT", "TCGAA", None, "ATCCG", None),
        ("GCTT", "GTTAC", "CAGGT", "TCGAG", "TACGG"),
        ("GCTT", "GTTAC", "CAGGT", "TCGAGTAA", "TACGG"),  # Test for stop codon
    ],
)
def test_make_construct(five_prime, three_prime, leader, c_gene, linker):
    cfg, airr_data, full = create_dummy_construct(
        five_prime, three_prime, leader, c_gene, linker
    )
    res = tcr_primers.assemble_one_construct(airr_data, ("VJ_1", "VDJ_1"), cfg)
    cons = res["sequence"]
    assert len(cons) == len(full)
    assert cons == full


@pytest.mark.parametrize(
    "five_prime,three_prime,leader,c_gene,linker",
    [
        (None, None, None, None, None),
        ("ATCGG", None, "GCTAT", "TAGCA", None),
        (None, "CGATC", "TTGAC", None, "ACGTA"),
        ("GGCAT", "TCGAA", None, "ATCCG", None),
        ("GCTT", "GTTAC", "CAGGT", "TCGAG", "TACGG"),
        ("GCTT", "GTTAC", "CAGGT", "TCGAGTAA", "TACGG"),  # Test for stop codon
    ],
)
def test_make_construct_class(five_prime, three_prime, leader, c_gene, linker):
    cfg, airr_data, full = create_dummy_construct(
        five_prime, three_prime, leader, c_gene, linker
    )
    c = tcr_cons.TCRConstruct(airr_data, ("VJ_1", "VDJ_1"), cfg)
    res = c.assemble()
    cons = res["sequence"]
    assert len(cons) == len(full)
    assert cons == full


# # TODO: finish this...
# @pytest.mark.parametrize("dct", [{}, {}])
# def test_validate_construct(dct):
#     tcr_primers.validate_construct()
