#!/usr/bin/env bash

# Helper script to clean up and prepare vcf file for standardization
clear_flag="${1}"
bcftools annotate -x "${clear_flag}" "${2}" | \
    awk 'BEGIN { FS="\t"; OFS="\t" } $5 !~ /\*/ {print}' | \
bcftools view -O z > tmp.vcf.gz
