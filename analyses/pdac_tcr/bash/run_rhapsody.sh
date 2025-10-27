#!/usr/bin/env bash

# In the rhapsodyPipeline directory
out="/data/project/stemcell/shannc/output/PDAC_TCR_2025-09-24/TCR1"
config="/data/home/shannc/chula-stem/analyses/pdac_tcr/rhapsody_1.yml"
# srun --qos=cpu24h --mem=50G \
#     --cpus-per-task=8 \
#     ./rhapsody pipeline --outdir "${out}" "${config}"

# [2025-10-17 Fri] try to count umis and cell indices
# See https://bd-rhapsody-bioinfo-docs.genomics.bd.com/resources/extra_utilities.html#annotate-cell-label-and-umi-only
# Ah forget about it. You need docker...
config="/data/home/shannc/chula-stem/analyses/pdac_tcr/bash/rhapsody_align_1.yml"
cwl-runner cwl/AnnotateCellLabelUMI.cwl "${config}"

# out="/data/project/stemcell/shannc/output/PDAC_TCR_2025-09-24/TCR2"
# config="/data/home/shannc/chula-stem/analyses/pdac_tcr/rhapsody_2.yml"

# srun --qos=cpu24h --mem=50G \
#     ./rhapsody pipeline --outdir "${out}" "${config}"
