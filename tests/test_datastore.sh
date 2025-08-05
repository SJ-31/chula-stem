#!/usr/bin/env bash

ds_actions -a record_cases -r ~/Bio_SDD/stem_synology/cancer_ngs/

./datastore/dummy_data.sh scrna_seq1 CN

datastore_actions -a store_away -r ~/Bio_SDD/stem_synology/cancer_ngs/ \
    -s ./datastore/new_cond.yaml \
    -o ./datastore/new_cond.csv

./datastore/dummy_data.sh scrna_seq2 BB

ds_actions -a store_away -r ~/Bio_SDD/stem_synology/cancer_ngs/ \
    -s ./datastore/existing_cond.yaml \
    -o ./datastore/existing_cond.csv
