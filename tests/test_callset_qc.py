import pytest
from chula_stem.callset_qc import _qc_main, merge_variant_calls


@pytest.mark.skip(reason="Done")
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
        min_VAF=0.05,
        accepted_filters="PASS",
        impact=False,
        canonical=True,
        informative=True,
        min_callers=1,
        tool_source_tag="SOURCE",
        vaf_adaptive=True,
    )
