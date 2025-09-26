#!/usr/bin/env bash

# In the rhapsodyPipeline directory
out="/data/project/stemcell/shannc/output/PDAC_TCR_2025-09-24/TCR1"
config="/data/home/shannc/chula-stem/analyses/pdac_tcr/rhapsody_1.yml"
srun --qos=cpu24h --mem=50G \
    ./rhapsody pipeline --outdir "${out}" "${config}"

out="/data/project/stemcell/shannc/output/PDAC_TCR_2025-09-24/TCR2"
config="/data/home/shannc/chula-stem/analyses/pdac_tcr/rhapsody_2.yml"

srun --qos=cpu24h --mem=50G \
    ./rhapsody pipeline --outdir "${out}" "${config}"
