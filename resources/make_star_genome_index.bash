#!/usr/bin/env bash

# Build star index for alignment in rnaseq
# Index generated <2025-01-09 Thu>
outdir="/data/project/stemcell/shannc/reference/tool_specific/GRCh38_star_index"
genome_fasta="/data/project/stemcell/shannc/reference/genomes/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
genome_gtf="/data/project/stemcell/shannc/reference/genomes/Homo_sapiens.GRCh38.113.gtf"
read_length=151

STAR \
    --runMode genomeGenerate \
    --genomeDir "${outdir}" \
    --genomeFastaFiles "${genome_fasta}" \
    --sjdbGTFfile "${genome_gtf}" \
    --sjdbOverhang "$((read_length - 1))"
