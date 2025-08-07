#!/usr/bin/env bash

config=$(realpath ./msi_config.yaml)
cd ../../snakemake/wes_msi_plot/ || exit
snakemake -c 1 --configfile "${config}"
