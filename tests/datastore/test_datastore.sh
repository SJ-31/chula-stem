#!/usr/bin/env bash

ds_actions -a record_cases -r ~/Bio_SDD/stem_synology/cancer_ngs/

./dummy_data.sh scrna_seq1 CN

ds_actions -a store_away -r ~/Bio_SDD/stem_synology/cancer_ngs/ \
    -s ./new_cond.yaml \
    -o ./new_cond.csv \
    -i

# ./dummy_data.sh scrna_seq2 BB

# ds_actions -a store_away -r ~/Bio_SDD/stem_synology/cancer_ngs/ \
#     -s ./existing_cond.yaml \
#     -o ./existing_cond.csv \
#     -n
