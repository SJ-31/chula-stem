#!/usr/bin/env bash

#SBATCH --job-name=format_variants
#SBATCH --cpus-per-task=1
#SBATCH --qos=cpu24h
#SBATCH --mem=30G

# Ran on Wed Dec  4 09:25:27 2024
# Format variant files for use with the WES and WES tumor-only pipelines
rename=$(yq ".rename.remove_chr" meta.yaml)
dir=$(yq ".dir.ref" meta.yaml)
dbsnp_renamed="${dir}/variants/dbSNP_renamed.vcf.gz"
gnomad_dir="${dir}/variants/gnomADv4.1.0_Exomes"

# Filter dbSNP to retain only known germline variants
dbsnp_g="${dir}/variants/dbSNP_renamed_germline.vcf.gz"
if [[ ! -e "${dbsnp_g}" ]]; then
    bcftools filter -i "INFO/SAO != 2" "${dbsnp_renamed}" -O z > "${dbsnp_g}"
    bcftools index "${dbsnp_g}"
fi

# Combine and rename gnomad files
all_gnomad="${dir}/variants/gnomADv4.1.0_all.vcf.gz"
if [[ ! -e "${all_gnomad}" ]]; then
    bcftools concat "${gnomad_dir}"/*.vcf.bgz | \
        bcftools annotate --rename-chrs "${rename}" -O z > "${all_gnomad}"
    bcftools index "${all_gnomad}"
fi

# Create GATK pileup file
pileup="${dir}/variants/gnomADv4.1.0_all_biallelic.vcf.gz"
if [[ ! -e "${pileup}" ]]; then
    ./prepare_pileup.bash "${all_gnomad}" "${pileup}"
    bcftools index "${pileup}"
fi
