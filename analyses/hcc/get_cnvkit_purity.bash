#!/usr/bin/env bash

dir="/data/project/stemcell/shannc/output/HCC/Exome/P17/variant_calling/5-P17-Cnvkit"
pushd "${dir}" || exit 1
vcf="/data/project/stemcell/shannc/output/HCC/Exome/P17/annotations/7-P17-Small_high_conf.vcf.gz"

for i in $(seq 1 9); do
    out="0.${i}_purity_call.cns"
    cnvkit.py call "5-P17-Cnvkit.call.cns" --vcf "${vcf}" --method clonal \
        --purity "0.${i}" \
        --output "${out}"
done
popd
