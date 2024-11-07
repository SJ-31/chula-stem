#!/usr/bin/env ipython

import vcfpy
import click

# Goal: make a cli command that adds an arbitrary tag TAG to vcf file with default value
# to every entry in the file
tag_name = "SOURCE"
description = "foo bar"
number = "."
type = "String"


def vcf_add_tag(name: str, description: str, number: str, type: str, default: str):
    new_line = vcfpy.OrderedDict(
        [("ID", name), ("Description", description), ("Number", number), ("Type", type)]
    )


sample = "/home/shannc/Bio_SDD/chula-stem/tests/sample.vcf.gz"
new_line = vcfpy.OrderedDict(
    [("ID", tag_name), ("Description", description), ("Number", number), ("Type", type)])
reader = vcfpy.Reader.from_path(sample)
# reader.add_info_line()
