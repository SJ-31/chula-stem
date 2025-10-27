#!/usr/bin/env bash

# TODO: [2025-10-04 Sat] waiting on the mixcr rep to help you
mixcr analyze bd-sc-xcr-rhapsody-full-length \
    --species hsa \
    --sample-sheet samples1.tsv \
    --tag-pattern "(R1:*)\(R2:*)\(INDEX1:)" \
    /data/project/stemcell/shannc/PDAC/raw/TCR_2025-09-23/TCRindexPCRno1_1.fastq.gz \
    /data/project/stemcell/shannc/PDAC/raw/TCR_2025-09-23/TCRindexPCRno1_2.fastq.gz \
    /data/project/stemcell/shannc/PDAC/raw/TCR_2025-09-23/SampleTAGindex1_1.fastq.gz \
    /data/project/stemcell/shannc/PDAC/raw/TCR_2025-09-23/SampleTAGindex1_2.fastq.gz \
    /data/project/stemcell/shannc/output/PDAC_TCR_2025-09-24/TCR1_MiXCR
