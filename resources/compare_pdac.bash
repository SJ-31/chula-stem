#!/usr/bin/env bash

# <2025-01-02 Thu> Script for comparing pdac runs
#
# Thu Jan  2 12:14:23 2025 Script is ready, just need to wait for pdac to finish to run
outdir="/data/project/stemcell/shannc/output/PDAC_comparison"
newdir="/data/project/stemcell/shannc/output/PDAC"
olddir="/data/project/stemcell/PDAC/Exome_PDAC_HN00215149"

pushd "${olddir}"
readarray -t patient_ids < <(ls PDAC* -d | sed 's/PDAC//')
popd

for id in ${patient_ids[@]}; do
    workdir="${outdir}/PDAC${id}"
    old_renamed="${workdir}/PDAC${id}.final.vcf.gz"
    new="${newdir}/P${id}/annotations/7-P${id}-Small_high_conf.vcf.gz"
    evaldir="${workdir}/vcfeval"
    new_moved="${workdir}/PDAC${id}_Small_high_conf.vcf.gz"

    if [[ ! -e "${new}" ]]; then
        continue
    fi

    if [[ ! -e "${new_moved}" ]]; then
        cp "${new}" "${new_moved}"
        tabix "${new_moved}"
    fi

    if [[ ! -e "${old_renamed}" ]]; then
        mkdir "${workdir}"
        old="${olddir}/PDAC${id}/PDAC${id}.final.vcf"
        ~/chula-stem/nextflow/bin/remove_chr.bash "${old}" "${old_renamed}"
        tabix "${old_renamed}"
    fi

    if [[ ! -e "${evaldir}" ]]; then
        rtg vcfeval \
        --bed-regions=/data/project/stemcell/shannc/reference/exome_kits/SureSelectHumanAllExonV6Hg38/Regions.bed.gz  \
        --output="${evaldir}" \
        --baseline="${old_renamed}" \
        --output-mode=combine \
        --calls="${new_moved}" \
        --template=/data/project/stemcell/shannc/reference/genomes/GRCh38.p14_filtered.sdf
    fi

done
