#!/usr/bin/env bash
dir="/data/project/stemcell/shannc/output/HCC/Exome/P17/variant_calling/5-P17-Cnvkit"
vcf="/data/project/stemcell/shannc/output/HCC/Exome/P17/annotations/7-P17-Small_high_conf.vcf.gz"
mutect_stats="/data/project/stemcell/shannc/output/HCC/Exome/P17/variant_calling/5-Mutect2/5-P17-Mutect2.vcf.gz.stats"

datadir="/data/project/stemcell/shannc/output/HCC/Exome/P17/PureCN_ref"

cnr="/ssh:shannc@161.200.107.77:/data/project/stemcell/shannc/output/HCC/Exome/P17/variant_calling/2025-1-20_Cnvkit_new_ref/4-P17_tumor-recal.cnr"

genome="/data/project/stemcell/shannc/reference/genomes/GRCh38.p14_filtered.fasta"
mappability="/data/project/stemcell/shannc/reference/tool_specific/GCA_000001405.15_GRCh38_no_alt_analysis_set_100.bw"

pon="???" # <2025-01-20 Mon> Need to wait for the wes pon to finish

