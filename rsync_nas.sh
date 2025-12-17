#!/usr/bin/env bash

exome="/data/project/stemcell/shannc/output/cancer_ngs/exome"
rna_seq="/data/project/stemcell/shannc/output/cancer_ngs/rna_seq"
scrna_seq="/data/project/stemcell/shannc/output/cancer_ngs/scrna_seq"
sc_atac_seq="/data/project/stemcell/shannc/output/cancer_ngs/sc_atac_seq"

target_root="/volume1/cancer_ngs"
port="22"

# rsync --archive --verbose -e "ssh -p ${port}" \
#     "${exome}" "stem_cell_admin@10.128.196.216:${target_root}"

# Should be able to do sync everything by just
# [2025-12-13 Sat] Getting permission errors with the 'archive' flag
rsync --recursive --verbose -e "ssh -p ${port}" \
    "${exome}" "${rna_seq}" \
    "stem_cell_admin@10.128.196.216:${target_root}"
