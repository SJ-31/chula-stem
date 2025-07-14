#!/usr/bin/env bash

#SBATCH --job-name=format_variants
#SBATCH --cpus-per-task=1
#SBATCH --qos=cpu24h
#SBATCH --mem=30G

# [2025-06-18 Wed]
# Rename sequences in variant files to prefix with chr to use with alignments
# from previous analysis
rename=$(yq ".rename.add_chr" meta.yaml)
refdir=$(yq ".dir.ref" meta.yaml)

gnomad_subset="$refdir/variants/gnomADv4.1.0_Exomes/random/gnomADv4.1_subset.vcf.gz"
gnomad_subset_biallelic="$refdir/variants/gnomADv4.1.0_Exomes/biallelic/all.vcf.gz"
dbsnp="$refdir/variants/dbSNP_renamed_germline.vcf.gz"

rename () {
    echo "${1}"
    bcftools annotate --rename-chrs "${rename}" -O z "${1}" > "${2}"
    if [[ -e "${2}" ]]; then
        gatk IndexFeatureFile -I "${2}"
        tabix "${2}"
    else
        echo "Bcftools failed to produce ${2}!"
        exit 1
    fi
}

gnomad_subset_chr="$refdir/variants/gnomADv4.1.0_Exomes/random/gnomADv4.1_subset_chr.vcf.gz"
gnomad_subset_biallelic_chr="$refdir/variants/gnomADv4.1.0_Exomes/biallelic/all_chr.vcf.gz"
dbsnp_chr="$refdir/variants/dbSNP_renamed_germline_chr.vcf.gz"

rename "${gnomad_subset}" "${gnomad_subset_chr}"
rename "${gnomad_subset_biallelic}" "${gnomad_subset_biallelic_chr}"
rename "${dbsnp}" "${dbsnp_chr}"
