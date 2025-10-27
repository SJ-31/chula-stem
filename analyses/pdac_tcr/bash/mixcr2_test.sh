#!/usr/bin/env bash

srun --qos=cpu24h --mem=50G \
    mixcr -Xmx40g analyze bd-sc-xcr-rhapsody-full-length-enhanced-bead-v2 \
    --species hsa \
    --sample-sheet ./mixcr_samplesheet_test.tsv \
    /data/project/stemcell/shannc/PDAC/raw/TCR_2025-09-23/TCRindexPCRno2_1.fastq.gz \
    /data/project/stemcell/shannc/PDAC/raw/TCR_2025-09-23/TCRindexPCRno2_2.fastq.gz \
    /data/project/stemcell/shannc/output/PDAC_TCR_2025-09-24/MiXCR_TCR2_test_v3/TCR2
