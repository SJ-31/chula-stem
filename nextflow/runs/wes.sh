#!/usr/bin/env bash

nextflow main.nf \
    --input ./input/test_crc_wes.csv \
    --output /data/project/stemcell/shannc/output \
    --reference_file ./input/wes_reference.nf \
    --routine "wes"
