#!/usr/bin/env bash

config=$(realpath ./env.yaml)
cd ../../snakemake/tcr/ || exit
report=../../analyses/output/pdac_tcr/report.zip
snakemake -c 1 --configfile "${config}" --report-after-run \
    --report "${report}" "$@"
unzip -q -o "${report}" -d ../../analyses/output/pdac_tcr
cd - || exit
