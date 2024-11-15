#!/usr/bin/env bash

#SBATCH --job-name=prepare_genome
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20G

# Prepares a genome fasta file to be used with the pipeline
genome="$1"
basename="${genome//.fasta/}"
gatk CreateSequenceDictionary -R "$genome"
bwa-mem2 index "$genome"
samtools faidx "$genome"

msi="msihomopoly-${basename}.tsv"
msisensor-pro scan -d "$genome" -o "$msi"
