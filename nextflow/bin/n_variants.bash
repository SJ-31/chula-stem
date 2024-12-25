#!/usr/bin/env bash

# Count the number of variants in input vcf file, accounting for potential duplicates
bcftools query -f "%CHROM %POS %ALT %REF" "${1}" | sort | uniq | wc -l
