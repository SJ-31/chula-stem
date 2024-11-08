#!/usr/bin/env ipython
# import chula_stem as

import vcfpy
from chula_stem.utils import _vcf_info_add_tag
from click.testing import CliRunner

runner = CliRunner()


def test_add_tag():
    test_dir = "/home/shannc/Bio_SDD/chula-stem/tests"
    sample = f"{test_dir}/sample.vcf.gz"
    output = f"{test_dir}/added.vcf.gz"
    _vcf_info_add_tag(
        "SOURCE", "Variant caller", ".", "String", "strelka2", sample, output
    )
    reader = vcfpy.Reader.from_path(output)
    for r in reader:
        r: vcfpy.Record
        try:
            assert r.INFO["SOURCE"][0] == "strelka2"
        except KeyError:
            continue
