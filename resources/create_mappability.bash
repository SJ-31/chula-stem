#!/usr/bin/env bash

# Creates a mappability map from a genome fasta file for use with delly
# Based on instructions by https://github.com/dellytools/delly
# WARNING: this might take a REALLY long time
#SBATCH --job-name=mappability_map
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --qos=cpu24h
genome="$1"
dicey chop --fq1 R1 --fq2 R2 "$genome" # Chops reference into PE reads
bwa-mem2 index "$genome"
bwa-mem2 mem R1.fq.gz R2.fq.gz | samtools sort -@ 8 -o srt.bam -
samtools index srt.bam
dicey mappability2 srt.bam # Compute mappability
gunzip map.fa.gz && bgzip map.fa && samtools faidx map.fa.gz
