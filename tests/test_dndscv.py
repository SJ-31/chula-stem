from chula_stem.dndscv import _dndscv


def test_dndscv():
    table = "/data/home/shannc/chula-stem/tests/dndscv_mutants.tsv"
    refdb = "/data/project/stemcell/shannc/reference/dndscv_GRCh38.rda"
    _dndscv(table, refdb, sample_id="sample2")
