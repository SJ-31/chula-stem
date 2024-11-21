#!/usr/bin/env bash


# Filter NCBI human genome to remove all but canonical chromosomes, then rename as only the chromosome number
input="$1"
output="$2"
cat "$input" | \
    seqkit grep -n -v -r \
        -p "alternate locus" \
        -p "unlocalized" \
        -p "mitochondrion" \
        -p "unplaced" \
        -p "genomic patch" | \
    sed -E 's/>.*chromosome ([0-9XY]+),.*/>\1/' > $output
