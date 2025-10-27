#!/usr/bin/env bash

config=$(realpath ./env.yaml)
cd ../../snakemake/tcr/ || exit
snakemake -c 1 --configfile "${config}" "$@"
cd - || exit
