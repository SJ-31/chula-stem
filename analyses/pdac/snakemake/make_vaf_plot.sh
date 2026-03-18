#!/usr/bin/env bash

date="${1}"

if [[ -z "${date}" ]]; then
    date=no_date
fi

config=$(realpath ./vaf_config.yaml)
cd ../../../snakemake/wes_vaf_plot/ || exit
snakemake -c 1 --configfile "${config}" --config date="${date}"
cd - || exit
