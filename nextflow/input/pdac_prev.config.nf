// * Required params
// Routine for processing the *recal files produced previously
workDir = "/data/project/stemcell/shannc/output/nf-work/PDAC"
includeConfig "$projectDir/input/wes_tumor_only_tools.config.nf"
params.outdir = "/data/project/stemcell/shannc/output/PDAC"
if (params.debug) {
    params.input = "$projectDir/input/manifests/pdac_DEBUG.csv"
} else {
    params.input = "$projectDir/input/manifests/pdac_prev.csv"
}
params.logdir = "/data/project/stemcell/shannc/output/PDAC/log"
params.routine = "wes"
params.cohort = "cohort"
params.tumor_only = true

// Map for pipeline resources
// See routine-specific documentation to see which need to be included

includeConfig "${projectDir}/input/references.config.nf"
params.report_data = ["hospital": "Dr.Nipan CU",
                      "physician": "Assoc. Prof. Nipan Israsena, M.D., Ph.D.",
                      "diagnosis": "Pancreatic ductal adenocarcinoma",
                      "sample_type": "organoid",
                      "collection_date": null,
]
params.report_text = ["front_page_details": "${projectDir}/reports/wes_front_page_details.txt",
                      "title": "Whole-exome sequencing results",
                      "disclaimer": "${projectDir}/reports/wes_disclaimer.txt",
                      "front_page_note": "${projectDir}/reports/wes_front_page_note.txt"]

genome_chr = "${refdir}/prev_assembly/" // TODO: wait for Aj to send it
targets_chr = "${kitdir}/S07604514_Regions_1.bed.gz"
baits_chr = "${kitdir}/S07604514_Covered_1.bed.gz"
baits_chr_unzipped = "${kitdir}/S07604514_Covered.bed"
genome_gff_chr = "${refdir}/prev_assembly/hg38.ensGene.gtf"

// TODO: gatk BedToIntervalList needs sequence dictionary, which you need to make
// from genome_chr
targets_il_chr = "${kitdir}/TODO"
baits_il_chr = "${kitdir}/TODO"

// TODO: write the file to make this in ~/Bio_SDD/chula-stem/resources/wes_format_variants_chr.bash
dbsnp_chr = "${refdir}/"//
gnomad_subset_chr = "${refdir}"
gnomad_subset_biallelic_chr = "${refdir}"
homopolymers_microsatellites_chr = "${refdir}"

params.ref = ["genome": genome_chr,
              "homopolymers_microsatellites": homopolymers_microsatellites,
              "targets": targets_chr,
              "baits": baits_chr,
              "baits_il": baits_il_chr,
              "targets_il": targets_il_chr,
              "baits_unzipped": baits_chr_unzipped,
              "genome_blacklist": genome_blacklist,
              "known_variants": [dbsnp_chr, gnomad_subset_chr],
              "genome_gff": genome_gff_chr,
              "germline": gnomad_subset_biallelic_chr,
              "dbsnp": dbsnp_chr,
              "pileup": gnomad_subset_biallelic_chr,
              "cnv_reference": cnv_reference,
              "msi_reference": msi_reference,
              "clingen_dosage": clingen_dosage,
              "civic_cache": null,
              "pandrugs2_cache": null,
              "cosmic_reference": cosmic_reference,
              "clingen_gene": clingen_gene]

// * Variant calling parameters

// ** Quality control filters
params.small_qc = ["accepted_filters": ["PASS"],
                   "min_tumor_depth": 10,
                   "max_normal_depth": 10,
                   "min_vaf": 0.10]
params.sv_qc = ["accepted_filters": ["PASS"]]

// * TEMPORARY
// <2025-01-08 Wed> Having difficulties with the therapy db caching, so for now
// will produce report without that data
process {
    withName: "REPORT" {
        ext.misc = ["plot": false, "show_therapies": true]
        ext.extra_paths = ["therapy_data_cache": "/data/home/shannc/.cache/2025-01-13_therapy_data.parquet"]
    }
}

// * Misc
// Add gnomad subsets to known_variants
gnomad_subsets = [
    "$refdir/variants/gnomADv4.1.0_Exomes/random/chr10.vcf.bgz",
    "$refdir/variants/gnomADv4.1.0_Exomes/random/chr11.vcf.bgz",
    "$refdir/variants/gnomADv4.1.0_Exomes/random/chr12.vcf.bgz",
    "$refdir/variants/gnomADv4.1.0_Exomes/random/chr13.vcf.bgz",
    "$refdir/variants/gnomADv4.1.0_Exomes/random/chr14.vcf.bgz",
    "$refdir/variants/gnomADv4.1.0_Exomes/random/chr15.vcf.bgz",
    "$refdir/variants/gnomADv4.1.0_Exomes/random/chr16.vcf.bgz",
    "$refdir/variants/gnomADv4.1.0_Exomes/random/chr17.vcf.bgz",
    "$refdir/variants/gnomADv4.1.0_Exomes/random/chr18.vcf.bgz",
    "$refdir/variants/gnomADv4.1.0_Exomes/random/chr19.vcf.bgz",
    "$refdir/variants/gnomADv4.1.0_Exomes/random/chr1.vcf.bgz",
    "$refdir/variants/gnomADv4.1.0_Exomes/random/chr20.vcf.bgz",
    "$refdir/variants/gnomADv4.1.0_Exomes/random/chr21.vcf.bgz",
    "$refdir/variants/gnomADv4.1.0_Exomes/random/chr22.vcf.bgz",
    "$refdir/variants/gnomADv4.1.0_Exomes/random/chr2.vcf.bgz",
    "$refdir/variants/gnomADv4.1.0_Exomes/random/chr3.vcf.bgz",
    "$refdir/variants/gnomADv4.1.0_Exomes/random/chr4.vcf.bgz",
    "$refdir/variants/gnomADv4.1.0_Exomes/random/chr5.vcf.bgz",
    "$refdir/variants/gnomADv4.1.0_Exomes/random/chr6.vcf.bgz",
    "$refdir/variants/gnomADv4.1.0_Exomes/random/chr7.vcf.bgz",
    "$refdir/variants/gnomADv4.1.0_Exomes/random/chr8.vcf.bgz",
    "$refdir/variants/gnomADv4.1.0_Exomes/random/chr9.vcf.bgz",
    "$refdir/variants/gnomADv4.1.0_Exomes/random/chrX.vcf.bgz",
    "$refdir/variants/gnomADv4.1.0_Exomes/random/chrY.vcf.bgz",
]

