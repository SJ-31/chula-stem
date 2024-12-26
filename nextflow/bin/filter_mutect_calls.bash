#!/usr/bin/env bash

while getopts "a:m:v:r:s:c:f:o:t:" f; do
    case "$f" in
        v) vcf=${OPTARG} ;;
        r) reference=${OPTARG} ;;
        t) segmentation=${OPTARG} ;;
        c) contamination=${OPTARG} ;;
        o) output=${OPTARG} ;;
        m) ro_model=${OPTARG} ;;
        s) stats=${OPTARG} ;;
        a) args=${OPTARG} ;;
        *)
            echo "Invalid flag"
            exit 1
            ;;
    esac
done


bcftools view "${vcf}" | \
    awk 'BEGIN { FS="\t"; OFS="\t" } $8 !~ /.*POPAF=[\.,a-zA-Z0-9]*inf;/ { print }' > tmp.vcf

cp "${stats}" tmp.vcf.stats

gatk IndexFeatureFile -I tmp.vcf

if [[ -n "${args}" ]]; then
    gatk FilterMutectCalls -V tmp.vcf \
        "${args}" \
        --ob-priors "${ro_model}" \
        --tumor-segmentation "${segmentation}" \
        --contamination-table "${contamination}" \
        --reference "${reference}" \
        -O "${output}"
else
    gatk FilterMutectCalls -V tmp.vcf \
        --ob-priors "${ro_model}" \
        --tumor-segmentation "${segmentation}" \
        --contamination-table "${contamination}" \
        --reference "${reference}" \
        -O "${output}"
fi
