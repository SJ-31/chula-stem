#!/usr/bin/env bash

#SBATCH --qos=cpu24h
#SBATCH --mem=5G
#SBATCH -o "output-${d}.txt"
#SBATCH -e "errors-${d}.txt"
nextflow "$@"
