#!/usr/bin/env bash

# Script for aggregating the repeat and CNV data for cross-referencing purposes
# Format: a tsv file of chr, start, stop, gene_name, source, type, accession
# Accessions are from dbVar

get_from_somatic_sv () {
    awk 'BEGIN { FS = "\t"; OFS="\t" }
  {
    gene = gensub(/.*affected:(.*), Position.*/, "\\1", "g", $13)
    source = "dbVar:"$14
    print $1,$2,$3,gene,source,$11,$4
  }' "$1"
}

get_header () {
  echo -e "chr\tstart\tstop\tgenes\tsource\ttype\taccession" > "${1}"
}

svdir="$(yq ".dir.ref" meta.yaml)/variants/structural"
gnomad_sv="${svdir}/gnomad.v4.1.sv.sites.vcf.gz"
dbvar_clinical="${svdir}/nstd102_clinical_sv.csv"
somatic_sv="${svdir}/somatic_sv.bed"
final_msi="${svdir}/aggregated_msi.tsv"
final_cnv="${svdir}/aggregated_cnv.tsv"

# Repetitive elements data
clinvar_msi="${svdir}/clinvar_microsatellites_somatic.vcf.gz" #
dbvar_msi="${svdir}/dbvar_repetitive.bed"
# gnomad_repeat="${svdir}/gnomad.v4.1.repetitive.vcf.gz" It appears that there are no
# repetitive elements catalogued in "gnomad_sv"
get_header "$final_msi"
if [[ ! -e "${dbvar_msi}" ]]; then
    awk 'BEGIN {FS = "\t"; OFS="\t"}
      { if ($10 == "Curated" && $11 == "tandem duplication")
      { print } }' "${somatic_sv}" > "${dbvar_msi}"
fi
bcftools query -f "%CHROM\t%POS\tNA\t%INFO/GENEINFO\tClinVar\t%CLNVCSO\t%DBVARID" "${clinvar_msi}" | \
  awk 'BEGIN { FS="\t"; OFS="\t" } { $4 = gensub(/\|/, ",", "g", $4); print }' >> "$final_msi"
get_from_somatic_sv "${dbvar_msi}" >> "$final_msi"


# CNV data
gnomad_cnv="${svdir}/gnomad.v4.1.cnv.non_neuro.vcf.gz"
dbvar_cnvs="${svdir}/dbvar_cnvs_somatic.bed"
dbvar_clinical_cnvs="${svdir}/dbvar_clinical_cnvs.csv" # Don't actually need this
# because somatic.bed already contains nstd102 data
if [[ ! -e "${dbvar_cnvs}" ]]; then
    awk 'BEGIN {FS = "\t"; OFS="\t"}
        $11 ~ /copy number/ { print }' "${somatic_sv}" > "${dbvar_cnvs}"
fi
if [[ ! -e "${dbvar_clinical_cnvs}" ]]; then
    awk 'BEGIN {FS = ","; OFS="\t"}
          { if ($4 ~ /copy number/)
            { print $0 } }' "${dbvar_clinical}" > "${dbvar_clinical_cnvs}"
fi
get_header "$final_cnv"
bcftools query -f "%CHROM\t%POS\t%INFO/END\t%INFO/Genes\tgnomADv4.1\t%INFO/SVTYPE\tNA" "${gnomad_cnv}" >> "$final_cnv"
get_from_somatic_sv "${dbvar_cnvs}" >> "$final_cnv"
