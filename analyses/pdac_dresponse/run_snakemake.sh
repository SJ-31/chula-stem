#!/usr/bin/env bash

config=$(realpath ./env.yaml)
cd ../../snakemake/drug_sensitivity/ || exit
snakemake -c 1 --configfile "${config}" "$@"
cd - || exit
