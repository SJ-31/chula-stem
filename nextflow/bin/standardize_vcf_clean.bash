#!/usr/bin/env bash

# Helper script to clean up and prepare vcf file for standardization
while getopts "i:c:r:o:" f; do
    case "$f" in
    r) reference=${OPTARG} ;; # Reference genome fasta
    c) clear_flag=${OPTARG} ;;
    i) input=${OPTARG} ;;
    o) output=${OPTARG} ;;
    *)
        echo "Invalid flag!"
        exit 1
        ;;
    esac
done

clean_vcf() {
    awk 'BEGIN { FS="\t"; OFS="\t" } $5 !~ /\*/ {print}' |
        bcftools norm --fasta-ref "${reference}" \
            --check-ref w -O z >"${1}"
}

if [[ -n "${clear_flag}" ]]; then
    bcftools annotate -x "${clear_flag}" "${input}" | clean_vcf "${output}"
else
    bcftools view "${input}" | clean_vcf "${output}"
fi
