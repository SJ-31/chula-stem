#!/usr/bin/env bash

dir="/data/project/stemcell/shannc/output/HCC/Exome/P17/variant_calling/5-P17-Cnvkit"
vcf="/data/project/stemcell/shannc/output/HCC/Exome/P17/annotations/7-P17-Small_high_conf.vcf.gz"
# The germline sample needs to contain both somatic and germline variants, so you need
# a pon

mutect_stats="/data/project/stemcell/shannc/output/HCC/Exome/P17/variant_calling/5-Mutect2/5-P17-Mutect2.vcf.gz.stats"
bam="/data/project/stemcell/shannc/output/HCC/Exome/P17/tumor/4-P17_tumor-recal.bam"
b1="/data/project/stemcell/shannc/output/WES_PON/HCC_19/normal/4-HCC_19_normal-recal.bam"
b2="/data/project/stemcell/shannc/output/WES_PON/HCC_8/normal/4-HCC_8_normal-recal.bam"

normal_outdir="/data/project/stemcell/shannc/output/HCC/Exome/P17/PureCN_normals"

datadir="/data/project/stemcell/shannc/output/HCC/Exome/P17/PureCN_ref"
datadir2="/data/project/stemcell/shannc/output/HCC/Exome/P17/PureCN_ref_cnvkit"
baits="/data/project/stemcell/shannc/reference/exome_kits/SureSelectHumanAllExonV6Hg38/Unzipped_covered.bed"
genome="/data/project/stemcell/shannc/reference/genomes/GRCh38.p14_filtered.fasta"
mappability="/data/project/stemcell/shannc/reference/tool_specific/GCA_000001405.15_GRCh38_no_alt_analysis_set_100.bw"

pon="/data/project/stemcell/shannc/output/WES_PON/tmp_pon.vcf.gz"

# 1. Get bait intervals
opt_baits="${datadir}/baits_optimized_hg38.bed"
bait_intervals="${datadir}/baits_intervals.txt"
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

cov_normalized="${datadir}/gc_normalized_cov.tsv"
cov_output="${datadir2}/4-P17_tumor-recal_coverage.txt.gz"
if [[ ! -e "${cov_output}" ]]; then
    Rscript "${PURECN}/Coverage.R" --out-dir "${datadir2}" \
        --bam "${bam}" \
        --intervals "${bait_intervals}"
fi

normal_bams=("${b1}" "${b2}")
for bm in "${normal_bams[@]}"; do
    Rscript "${PURECN}/Coverage.R" --out-dir "${normal_outdir}" \
        --bam "${bm}" \
        --intervals "${bait_intervals}"
done

# <2025-02-27 Thu> Need NORMAL coverages
coverages="${datadir2}/coverage_list.txt"
ls ${normal_outdir}/*_coverage.txt.gz > "${coverages}"

# 3. Build normal db (this requires the coverage files from normals generated above)
# coverages=
# This creates the normal db AND the mapping bias file, both of which are RDS
Rscript "${PURECN}/NormalDB.R" --out-dir ${datadir} \
    --coverage-files "${coverages}" \
    --normal-panel "${pon}" \
    --genome hg38 \
    --assay agilent_v6

# 4. Run PureCN to normalize, segment and get purity + ploidy
# Rscript "${PURECN}/PureCN.R" --out "${data}" \
    #     --tumor "${}" \ # Tumor coverage file
#     --sampleid P17 \
#     --vcf "${vcf}" \ # tumor vcf
#     --stats-file "${mutect_stats}" \
#     --normaldb "${}" \
#     --intervals "${bait_intervals}" \
#     --genome hg38
#
