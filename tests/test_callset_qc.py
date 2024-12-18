import pytest
from chula_stem.callset_qc import merge_variant_calls


def test_variant_merge():
    input = (
        "/home/shannc/Bio_SDD/chula-stem/tests/vep/7-patient_10-VEP_small_1.tsv"
    )
    result = merge_variant_calls(
        input,
        tool_source_tag="SOURCE",
        minimum_callers=1,
        grouping_cols=["Loc", "Ref", "Alt"],
    )
    result.write_csv("test_variant_merge1.tsv", separator="\t")


def _test_main():
    input = (
        "/home/shannc/Bio_SDD/chula-stem/tests/vep/7-patient_10-VEP_small_1.tsv"
    )
