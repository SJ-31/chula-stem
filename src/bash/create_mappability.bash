#!/usr/bin/env bash

# Instructions taken from PureCN vignette at https://www.bioconductor.org/packages/release/bioc/vignettes/PureCN/inst/doc/PureCN.pdf

reference="$1" # Genome fasta file
kmer="$2" # Kmer size should be the length of mapped reads
pref=${reference/\.*/}
threads=4
gem-indexer -T ${threads} -c dna -i "${reference}" -o "${pref}"_index
gem-mappability -T ${threads} -i "${pref}_index.gem" -l "${kmer}" \
-o "${pref}_${kmer}" -m 2 -e 2
gem-2-wig -i "${pref}_index.gem" -i "${pref}_${kmer}.mappability" \
