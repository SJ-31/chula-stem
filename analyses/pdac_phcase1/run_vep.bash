#!/usr/bin/env bash
input="./725-69_PDAC_organoid_PHcase1_DNA_10GB_725-71_PHcase1_from_PBMC_DNA_10GB.mutect2.vcf.gz"

singularity exec ~/tools/vep.sif vep --cache \
    --input_file ${input} \
    --clin_sig_allele 0 \
    --pubmed \
    --overlaps \
    --exclude_predicted \
    --gene_phenotype \
    --regulatory \
    --canonical \
    --offline \
    --dir_cache /data/home/shannc/.cache/vep \
    --species homo_sapiens \
    --stats_file vep_stats.html \
    --hgvs \
    --hgvsg \
    --vcf \
    --vcf_info_field ANN \
    --compress_output bgzip \
    --output_file PHase1_annotated.vcf.gz
