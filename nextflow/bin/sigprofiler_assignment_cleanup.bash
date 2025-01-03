#!/usr/bin/env bash

# 1. Remove ambiguous alt alleles "*"
# 2. Remove duplicate variant calls (likely due to multiple tools)
renamed=$(basename "${1}" .gz)
bcftools view "${1}" | \
    awk 'BEGIN { FS="\t"; OFS="\t" } $5 !~ /\*/ {print}' | \
    bcftools norm --rm-dup exact > "${2}/${renamed}"
