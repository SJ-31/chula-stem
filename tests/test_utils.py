#!/usr/bin/env ipython
# import chula_stem as

import vcfpy
from chula_stem.utils import vcf_info_add_tag


def test_add_tag():
    test_dir = "/home/shannc/Bio_SDD/chula-stem/tests"
    sample = f"{test_dir}/sample.vcf.gz"
    output = f"{test_dir}/added.vcf.gz"
    vcf_info_add_tag(
        "SOURCE", "Variant caller", ".", "String", "strelka2", sample, output
    )
    reader = vcfpy.Reader.from_path(output)
