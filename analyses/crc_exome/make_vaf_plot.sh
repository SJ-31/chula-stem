#!/usr/bin/env bash

config=$(realpath ./2026-05-22_vaf_config.yml)
cd ../../snakemake/wes_vaf_plot/ || exit
snakemake -c 1 --configfile "${config}" --keep-incomplete --rerun-triggers mtime
