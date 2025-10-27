#!/usr/bin/env bash

# [2025-10-21 Tue] Changing the default tag pattern based on the BD cls whitelist
# You added the V3-specific Tso_capture_seq to the end of the UMI, as described in the
# file https://bd-rhapsody-public.s3.amazonaws.com/CellLabel/rhapsody_cell_label.py.txt
# [2025-10-23 Thu] this didn't do anything
#
srun --qos=cpu24h --mem=50G --cpus-per-task=8 \
    mixcr -Xmx40g analyze bd-sc-xcr-rhapsody-full-length-enhanced-bead-v2 \
    --species hsa \
    --tag-pattern "^acaggaaactcatggtgcgt(CELL1:N{9})aatg(CELL2:N{9})ccac(CELL3:N{9})(UMI:N{8})gtggagtcgtgattata\^(R2:*)" \
    --threads 8 \
    --sample-sheet ./mixcr_samplesheet_test.tsv \
    /data/project/stemcell/shannc/PDAC/raw/TCR_2025-09-23/TCRindexPCRno2_1.fastq.gz \
    /data/project/stemcell/shannc/PDAC/raw/TCR_2025-09-23/TCRindexPCRno2_2.fastq.gz \
    /data/project/stemcell/shannc/output/PDAC_TCR_2025-09-24/MiXCR_TCR2_test_v3_with_tso/TCR2 "$@"
