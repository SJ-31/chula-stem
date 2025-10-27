#!/usr/bin/env bash

# [2025-10-21 Tue] Testing to see if using the imgt library (what BD uses) increases
# the number of called cells

srun --qos=cpu24h --mem=50G \
    mixcr -Xmx40g analyze bd-sc-xcr-rhapsody-full-length-enhanced-bead-v2 \
    --species hsa \
    --library imgt \
    --sample-sheet ./mixcr_samplesheet_test.tsv \
    /data/project/stemcell/shannc/PDAC/raw/TCR_2025-09-23/TCRindexPCRno2_1.fastq.gz \
    /data/project/stemcell/shannc/PDAC/raw/TCR_2025-09-23/TCRindexPCRno2_2.fastq.gz \
    /data/project/stemcell/shannc/output/PDAC_TCR_2025-09-24/MiXCR_TCR2_test_imgt/TCR2
