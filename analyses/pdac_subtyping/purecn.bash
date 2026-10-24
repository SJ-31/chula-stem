#!/usr/bin/env bash

# The germline sample needs to contain both somatic and germline variants, so you need
# a pon

dir="/data/project/stemcell/shannc/output"
refdir="/data/project/stemcell/shannc/reference"
outdir="/data/project/stemcell/shannc/tests/purecn_pdac_2026-03-06"

tbam="${dir}/PHCase/PHcase_1/tumor/4-PHcase_1_tumor-recal.bam"

normal_outdir="/data/project/stemcell/shannc/output/HCC/Exome/P17/PureCN_normals"

PURECN="/data/home/shannc/miniforge3/lib/R/library/PureCN/extdata"

baits="${refdir}/exome_kits/SureSelectHumanAllExonV6Hg38/Unzipped_covered.bed"
targets="${refdir}/exome_kits/SureSelectHumanAllExonV6Hg38/Regions.bed.gz"
genome="${refdir}/genomes/GRCh38.p14_filtered.fasta"
mappability="${refdir}/tool_specific/GCA_000001405.15_GRCh38_no_alt_analysis_set_100.bw"
gnomad_subset_biallelic="${refdir}/variants/gnomADv4.1.0_Exomes/biallelic/all.vcf.gz"

# Setup
mutect_out="${outdir}/mutect_with_germline.vcf.gz"

if [[ ! -e "${mutect_complete}" ]]; then
    gatk Mutect2 -R "${genome}" \
        -I "${tbam}" \
        --genotype-germline-sites true \
        --genotype-pon-sites true \
        --intervals "${targets}" \
        --germline-resource "${gnomad_subset_biallelic}" \
        --output "${mutect_out}" \
        --interval-padding 50
fi

# 1. Get bait intervals
opt_baits="${outdir}/baits_optimized_hg38.bed"
bait_intervals="${outdir}/baits_intervals.txt"
if [[ ! -e "${opt_baits}" ]]; then
    Rscript "${PURECN}/IntervalFile.R" --in-file "${baits}" \
        --fasta "${genome}" \
        --off-target --genome hg38 \
        --out-file "${bait_intervals}" \
        --export "${opt_baits}" \
        --mappability "${mappability}"
fi

# 2. Calculate coverage for each sample with PureCN
# if you have normals, use the same flags
# Rscript "${PURECN}/Coverage.R" --outdir "${datadir}" \
#     --bam "${bam}" \
#     --intervals "${bait_intervals}"

# Or use CNVKit's coverage file and let PureCN only GC-normalize
# <2025-01-20 Mon> this doesn't seem to be working correctly
# Need to give it the bam
# coverage="/data/project/stemcell/shannc/output/HCC/Exome/P17/variant_calling/2025-1-20_Cnvkit_new_ref/4-P17_tumor-recal.targetcoverage"

# cov_output="${outdir}/coverage.txt.gz"
# if [[ ! -e "${cov_output}" ]]; then
#     Rscript "${PURECN}/Coverage.R" --out-dir "${outdir}" \
#         --bam "${tbam}" \
#         --intervals "${bait_intervals}"
# fi
# The automatic name for this is FILENAME_coverage.txt.gz

# Rscript "${PURECN}/Coverage.R" --out-dir "${outdir}" \
#     --bam ./normals.list \
#     --intervals "${bait_intervals}"

normal_coverages="${outdir}/coverage_list.txt"
ls ${outdir}/*normal*_coverage.txt.gz >"${normal_coverages}"

# 3. Build normal db (this requires the coverage files from normals generated above)
# coverages=
# This creates the normal db AND the mapping bias file, both of which are RDS
# TODO: update this
Rscript "${PURECN}/NormalDB.R" --out-dir ${datadir} \
    --coverage-files "${coverages}" \
    --normal-panel "${pon}" \
    --genome hg38 \
    --assay agilent_v6

# 4. Run PureCN to normalize, segment and get purity + ploidy
# Rscript "${PURECN}/PureCN.R" --out . \
#     --tumor "${}" \ # Tumor coverage file
#     --sampleid P17 \
#     --vcf "${vcf}" \ # tumor vcf
#     --stats-file "${mutect_stats}" \
#     --normaldb "${}" \
#     --intervals "${bait_intervals}" \
#     --genome hg38
#
