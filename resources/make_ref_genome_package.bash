#!/usr/bin/env bash

# Create a reference genome package for use with accucopy

package_folder="$1"
ref="$2" # IMPORTANT: chromosomes of the reference must be labelled as chr1, chr2, ... not 1, 2
snps="$3" # The common SNP file in vcf format, will be converted to BED

b1=${snps/.vcf/}
basename=${b1/.gz/}
bed="$package_folder/${basename}.bed"

mkdir "$package_folder"
mv "$ref" "$package_folder"
vcf2bed "$snps" > "$bed"
samtools faidx "$package_folder/$ref"
gatk CreateSequenceDictionary -R "$package_folder/$ref"
tabix "$bed"
