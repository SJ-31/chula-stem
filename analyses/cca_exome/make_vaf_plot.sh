#!/usr/bin/env bash

config=$(realpath ./vaf_config.yaml)
cd ../../snakemake/wes_vaf_plot/ || exit
snakemake -c 1 --configfile "${config}"
