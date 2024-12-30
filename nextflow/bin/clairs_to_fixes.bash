#!/usr/bin/env bash

# Remove weird calls where there is no info field (and just a '.' instead)

bcftools view "${1}" | \
    awk 'BEGIN { FS="\t"; OFS="\t" }
   $8 !~ /^\..*/ { print }' | \
       bcftools view -O z > "${2}"
