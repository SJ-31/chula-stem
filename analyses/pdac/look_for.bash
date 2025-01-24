#!/usr/bin/env bash

dir="../data_all/output/PDAC"

for i in ${dir}/*/annotations/7-P*-VEP_small.vcf.gz; do
    count=$(bcftools query -f "[%AF %AD]\t%INFO/ANN" "${i}" | grep KRAS | cut -f 1)
    echo "SAMPLE $(basename "${i}")"
    echo "${count}"
    echo ""
done
