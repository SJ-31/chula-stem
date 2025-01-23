#!/usr/bin/env bash

genome="/data/project/stemcell/shannc/reference/genomes/GRCh38.p14_filtered.fasta"
bam="/data/project/stemcell/shannc/output/HCC/Exome/P17/tumor/4-P17_tumor-recal.bam"
baits="/data/project/stemcell/shannc/reference/exome_kits/SureSelectHumanAllExonV6Hg38/Unzipped_covered.bed"
access="/data/project/stemcell/shannc/output/HCC/Exome/cnvkit_cnn/4-flat_reference-CnvkitPrep_access.bed"
cnvkit.py autobin -m hybrid \
    -f "${genome}" \
    -t "${baits}" \
    -g "${access}" \
    "${bam}" > binned.txt
