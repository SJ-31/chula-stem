#!/usr/bin/env bash

# Script for retrieving only biallelic sites from a vcf file
# Intended for use with gatk GetPileupSummaries
#
vcf="$1"
output="$2"
if [[ -z $( bcftools head "${vcf}" | grep "ID=AF,") ]]; then
    echo "The given file has no population allele frequencies (AF) in the INFO field!"
    exit 1
fi
bcftools view -m2 -M2 -v snps "$1" > "$output"
# -m2 = min alleles 2
# -M2 = Max alleles 2
