#!/usr/bin/env bash

# [2025-10-23 Thu] You've observed that mixcr will drop all alignments that don't cover cdr3, of which there are many, unlike BD pipeline. Hopefully adding this flag will fix the issue

srun --qos=cpu24h --mem=50G --cpus-per-task=8 \
    mixcr -Xmx40g analyze bd-sc-xcr-rhapsody-full-length-enhanced-bead-v2 \
    --species hsa \
    --threads 8 \
    --keep-non-CDR3-alignments \
    --sample-sheet ./mixcr_samplesheet_test.tsv \
    /data/project/stemcell/shannc/PDAC/raw/TCR_2025-09-23/TCRindexPCRno2_1.fastq.gz \
    /data/project/stemcell/shannc/PDAC/raw/TCR_2025-09-23/TCRindexPCRno2_2.fastq.gz \
    /data/project/stemcell/shannc/output/PDAC_TCR_2025-09-24/MiXCR_TCR2_test_no_drop_cdr3/TCR2 "$@"
