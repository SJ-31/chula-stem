#!/usr/bin/env bash

# Split a fasta file into separate fasta files, one for each entry, storing the results in
# "outdir"

fasta="$1"
outdir="$2"
mkdir "${outdir}"
cd "${outdir}"
awk 'BEGIN{RS=">";FS="\n"} NR>1{fnme=$1".fasta"; print ">" $0 > fnme; close(fnme);}' "$fasta"
cd -
