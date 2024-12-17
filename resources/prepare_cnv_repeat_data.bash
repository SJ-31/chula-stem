#!/usr/bin/env bash

# Script for aggregating the repeat and CNV data for cross-referencing purposes
# The format of the final files is a tsv file of chr, start, stop, gene_name, source, type, clinical_status

svdir="$(yq ".dir.ref" meta.yaml)/variants/structural"
tdir="$(yq ".dir.ref" meta.yaml)/variants/therapeutics"
gnomad_sv="${svdir}/gnomad.v4.1.sv.sites.vcf.gz"
dbvar_clinical="${svdir}/nstd102_clinical_sv.csv"
somatic_sv="${svdir}/somatic_sv.bed"

# Repetitive elements data
clinvar_msi="${svdir}/clinvar_microsatellites_somatic.vcf.gz" #
dbvar_msi="${svdir}/dbvar_repetitive.bed"
# gnomad_repeat="${svdir}/gnomad.v4.1.repetitive.vcf.gz" It appears that there are no
# repetitive elements catalogued in "gnomad_sv"
if [[ ! -e "${dbvar_msi}" ]]; then
    awk 'BEGIN {FS = "\t"; OFS="\t"}
      { if ($10 == "Curated" && $11 == "tandem duplication")
      { print } }' "${somatic_sv}" > "${dbvar_msi}"
fi



# CNV data
# The format of the final files is a tsv file of chr, start, stop, gene_name, source, type, clinical_status
gnomad_cnv="${svdir}/gnomad.v4.1.cnv.non_neuro.vcf.gz"
dbvar_cnvs="${svdir}/dbvar_cnvs_somatic.bed"
dbvar_clinical_cnvs="${svdir}/dbvar_clinical_cnvs.csv"
if [[ ! -e "${dbvar_cnvs}" ]]; then
    awk 'BEGIN {FS = "\t"; OFS="\t"}
        $11 ~ /copy number/ { print }' "${somatic_sv}" > "${dbvar_cnvs}"
fi
if [[ ! -e "${dbvar_clinical_cnvs}" ]]; then
    awk 'BEGIN {FS = ","; OFS="\t"}
          { if ($4 ~ /copy number/)
            { print $0 } }' "${dbvar_clinical}" > "${dbvar_clinical_cnvs}"
fi
echo -e "chr\tstart\tstop\tgenes\tsource\ttype" > cnvs.tsv
bcftools query -f "%CHROM\t%POS\t%INFO/END\t%INFO/Genes\tgnomADv4.1\t%INFO/SVTYPE" "${gnomad_cnv}" > 1.tsv
awk '{}'


# TODO: aggregate the other three into a single file
