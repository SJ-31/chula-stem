#!/usr/bin/env bash

sample="$1"

bcftools annotate --rename-chrs ~/remote/reference/rename/remove_chr.tsv "7-${sample}-VEP_small_original.vcf.gz" |
    sed 's/TOOL_SOURCE/SOURCE/' >"7-${sample}-VEP_small.vcf.gz"
