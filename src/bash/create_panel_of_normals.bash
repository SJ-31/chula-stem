#!/usr/bin/env bash

while getopts "l:v:m:o:t:a:" f; do
    case "$f" in
        l) sample_vcfs=${OPTARG} ;; # Text file containing paths of in-house normal sample vcfs
        m) minimum_samples=${OPTARG} ;; # Variants only accepted if they are found in this
        # number of samples
        v) other_vcfs=${OPTARG};; # Text file listing paths of external vcfs  e.g. dbSNP to be added into the PON
        o) output_cleaned=${OPTARG} ;;
        a) output=${OPTARG} ;;
        t) threads=${OPTARG} ;;
        *) echo "Invalid flag"; exit 1;;
    esac
done

clean_vcf () {
    bcftools view -G | \
        bcftools annotate -x INFO,FILTER,FORMAT,QUAL -O z -W -o "${1}" --threads "${threads}"
}

mkdir cleaned

while read -r v; do
    bcftools view "${v}" | clean_vcf cleaned/"$(basename "${v}")"
done < "${other_vcfs}"

bcftools isec -l "${sample_vcfs}" -n "+${minimum_samples}" -O z -p tmp
bcftools merge tmp/*vcf.gz | clean_vcf cleaned/samples_merged.vcf.gz
bcftools concat cleaned/*.vcf.gz -D -a -O z > "${output_cleaned}"

bcftools merge tmp/*vcf.gz -O z > "${output}"

rm -R tmp
rm -R cleaned
