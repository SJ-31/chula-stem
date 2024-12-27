import pytest
from chula_stem.callset_qc import _qc_main, merge_variant_calls
from chula_stem.utils import contain_join
import polars as pl


@pytest.mark.skip(reason="Done")
def test_variant_merge():
    input = "/home/shannc/Bio_SDD/chula-stem/tests/vep/7-patient_10-VEP_small_1.tsv"
    result = merge_variant_calls(
        input,
        tool_source_tag="SOURCE",
        minimum_callers=1,
        grouping_cols=["Loc", "Ref", "Alt"],
    )
    result.write_csv("test_variant_merge1.tsv", separator="\t")


@pytest.mark.skip(reason="Done")
def test_main():
    input = "/home/shannc/Bio_SDD/chula-stem/tests/vep/format/6-null-VEP.tsv"
    # _qc_main(
    #     input,
    #     "/home/shannc/Bio_SDD/chula-stem/tests/vep/format/vep_qc.tsv",
    #     min_tumor_depth=10,
    #     max_normal_depth=5,
    #     min_VAF=0.05,
    #     accepted_filters="PASS",
    #     impact=True,
    #     canonical=True,
    #     informative=True,
    #     min_callers=1,
    #     tool_source_tag="SOURCE",
    #     vaf_adaptive=False,
    # )
    _qc_main(
        input,
        "/home/shannc/Bio_SDD/chula-stem/tests/vep/format/vep_qc2.tsv",
        min_tumor_depth=10,
        max_normal_depth=5,
        min_vaf=0.05,
        accepted_filters="PASS",
        impact=False,
        canonical=True,
        informative=True,
        min_callers=1,
        tool_source_tag="SOURCE",
        vaf_adaptive=True,
        ignore_regions="/home/shannc/Bio_SDD/chula-stem/tests/msisensor/4-null-CR.tsv",
        chr_col=1,
        start_col=9,
        end_col=10,
    )


def test_contain():
    d1 = pl.DataFrame(
        {
            "a": ["foo", "bar", "baz"],
            "b1": [21, 50, 50],
            "b2": [30, 70, 80],
            "d": ["a", "b", "c"],
        }
    )
    d2 = pl.DataFrame(
        {
            "a": ["foo", "bar", "baz"],
            "c1": [20, 50, 76],
            "c2": [40, 80, 90],
            "d": ["a", "f", "c"],
        }
    )
    ans1 = pl.DataFrame(
        {
            "a": ["foo", "bar"],
            "b1": [21, 50],
            "b2": [30, 70],
            "d": ["a", "b"],
            "c1": [20, 50],
            "c2": [40, 80],
            "d_right": ["a", "f"],
        }
    )
    ans2 = pl.DataFrame(
        {
            "a": ["foo"],
            "b1": [21],
            "b2": [30],
            "d": ["a"],
            "c1": [20],
            "c2": [40],
        }
    )
    result = contain_join(d1, d2, "b1", "c1", "c2", on=["a"], x_end="b2")
    result2 = contain_join(d1, d2, "b1", "c1", "c2", on=["a", "d"], x_end="b2")
    assert ans1.equals(result)
    assert ans2.equals(result2)
    print(result)
    print(result2)
