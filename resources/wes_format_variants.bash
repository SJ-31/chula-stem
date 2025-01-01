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
# Completed Fri Dec  6 08:39:35 2024
dbsnp_g="${dir}/variants/dbSNP_renamed_germline.vcf.gz"
if [[ ! -e "${dbsnp_g}" ]]; then
    echo "Starting germline filter"
    bcftools filter -i "INFO/SAO != 2" "${dbsnp_renamed}" -O z > "${dbsnp_g}"
    bcftools index "${dbsnp_g}"
else
    echo "dbSNP germline completed"
fi

# Randomly sample from gnomad variants and get biallelic
random="${dir}/variants/gnomADv4.1.0_Exomes/random"
biallelic="${dir}/variants/gnomADv4.1.0_Exomes/biallelic"
for g in "${dir}"/variants/gnomADv4.1.0_Exomes/*vcf.bgz; do
    basename=$(echo "${g}" | sed -r 's/.*(chr[0-9XY]+.*)/\1/')
    uncompressed=$(echo "${basename}" | sed 's/\.bgz//')

    if [[ ! -e "${random}/${basename}" ]]; then
       echo "Randomly sampling from ${g}"
       gatk IndexFeatureFile -I "${g}"
       gatk SelectVariants --variant "${g}" \
           --select-random-fraction 0.1 \
           --output "${random}/${uncompressed}"
       bcftools annotate --rename-chrs "${rename}" -O z "${random}/${uncompressed}" \
           > "${random}/${basename}"
       rm "${random}/${uncompressed}"
    else
        echo "${g} already sampled"
    fi

    if [[ ! -e "${biallelic}/${basename}" ]]; then
       echo "Getting biallelic sites from ${g}"
       gatk IndexFeatureFile -I "${random}/${basename}"
       gatk SelectVariants --variant "${random}/${basename}" \
           --output "${biallelic}/${uncompressed}" \
           --restrict-alleles-to BIALLELIC
       bcftools view -m2 -M2 -v snps "${random}/${basename}" -O z \
           > "${biallelic}/${basename}"
       gatk IndexFeatureFile -i "${biallelic}/${basename}"
    else
        echo "Biallelic from ${g} retrieved"
    fi
done

for dir in $("${biallelic}" "${random}"); do
    cd "${dir}"
    bcftools concat *.vcf.gz -O z > all.vcf.gz
    gatk IndexFeatureFile -i all.vcf.gz
done
