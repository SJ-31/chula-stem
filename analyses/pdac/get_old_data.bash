#!/usr/bin/env bash

dir="/data/project/stemcell/PDAC/Exome_PDAC_HN00215149"
outdir="/data/project/stemcell/shannc/output/PDAC/OLD"
tools="/data/home/shannc/tools"
pushd "${dir}"
for xls in */*.xlsx; do
    name=$(basename "${xls}" | sed 's/_SNP_Indel_ANNO.xlsx//')
    xlsx2csv "${xls}" "${outdir}/${name}.csv"
done
popd
