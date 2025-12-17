#!/usr/bin/env bash

public_data="/data/project/stemcell/shannc/public_data"
reference="/data/project/stemcell/shannc/reference"

target_root="/volume1/homes/shannc"
port="22"

rsync --recursive --verbose -e "ssh -p ${port}" \
    "${public_data}" "${reference}" \
    "stem_cell_admin@10.128.196.216:${target_root}"
