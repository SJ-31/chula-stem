#!/usr/bin/env ipython
# import chula_stem as
import pytest
from chula_stem.utils import _format_vep_vcf

testdir = "/home/shannc/Bio_SDD/chula-stem/tests/vep"


# @pytest.mark.skip(reason="Done")
def test_format_vep():
    sv = "/home/shannc/Bio_SDD/chula-stem/tests/vep/7-patient_10-VEP_SV.vcf.gz"
    small = (
        "/home/shannc/Bio_SDD/chula-stem/tests/vep/7-patient_10-VEP_small.vcf.gz"
    )
    _format_vep_vcf(
        sv,
        f"{testdir}/sv.tsv",
        tumor_sample="10_cancer",
        normal_sample="10_blood",
        tool_source_tag="SOURCE",
    )
    _format_vep_vcf(
        small,
        f"{testdir}/small.tsv",
        tumor_sample="10_cancer",
        normal_sample="10_blood",
        tool_source_tag="SOURCE",
    )
    _format_vep_vcf(
        small,
        f"{testdir}/small_no.tsv",
        tumor_sample="10_cancer",
        tool_source_tag="SOURCE",
    )
