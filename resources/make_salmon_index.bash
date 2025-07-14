#!/usr/bin/env bash
transcriptome="/data/project/stemcell/shannc/reference/transcriptomes/Homo_sapiens.GRCh38.cdna.all.fa"
genome="/data/project/stemcell/shannc/reference/genomes/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"

grep "^>" "${genome}" | cut -f 1 -d " " | sed 's/>//g' > decoys.txt
cat "${transcriptome}" "${genome}" > salmon_tmp.fa
salmon index -t salmon_tmp.fa -d decoys.txt -i salmon_index -p 8 -k 31
