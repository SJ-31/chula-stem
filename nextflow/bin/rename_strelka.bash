#!/usr/bin/env bash

while getopts "p:n:t:o:e:d:" f; do
    case "$f" in
        o) out=${OPTARG};;
        t) tumor=${OPTARG} ;;
        n) normal=${OPTARG} ;;
        e) source_name=${OPTARG} ;;
        d) source_description=${OPTARG} ;;
        p) prefix=${OPTARG} ;;
        *)
            echo "Flag not recognized"
            exit 1
            ;;
    esac
done

mv ${out}/results/variants/*.vcf.gz .
for variant in somatic*.vcf.gz; do
    base=$(echo "$variant" | sed -e 's/\.vcf\.gz//' -e 's/somatic\.//')
    name="${prefix}-${base}_Strelka.vcf"

    rename_vcf.bash -v -i $variant -n "${normal}" -t "${tumor}" | \
        vcf_info_add_tag.bash -n "${source_name}" \
            -d "${source_description}" \
            -b '.' \
            -t String \
            -a strelka2 \
            -o "${name}"

    rm "$variant"
done
