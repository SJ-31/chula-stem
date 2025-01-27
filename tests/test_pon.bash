#!/usr/bin/env bash

merged="/data/project/stemcell/shannc/tests/test_pon/merged.vcf.gz"
bcftools view "${merged}" | ../src/chula_stem/pon.py
