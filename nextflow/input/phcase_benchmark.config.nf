// * Required params
workDir = "/data/project/stemcell/shannc/output/nf-work/PHCase"
includeConfig "$projectDir/input/wes_tools.config.nf"
params.outdir = "/data/project/stemcell/shannc/output/Phcase_benchmark"
params.input = "$projectDir/input/manifests/pdac_phcase1.csv"
params.logdir = "${params.outdir}/log"
params.routine = "wes"
params.cohort = "PHcase"
params.tumor_only = false
params.source_name = "SOURCE"

report {
    enabled = true
    file = "${params.outdir}/report.html"
}

params.report_data = ["hospital": "Dr.Nipan CU",
                      "physician": "Assoc. Prof. Nipan Israsena, M.D., Ph.D.",
                      "diagnosis": "PDAC",
                      "sample_type": "organoid",
                      "collection_date": null,
]
params.report_text = ["front_page_details": "/data/home/shannc/chula-stem/nextflow/reports/wes_front_page_details.txt",
                      "title": "Whole-exome sequencing results",
                      "disclaimer": "/data/home/shannc/chula-stem/nextflow/reports/wes_disclaimer.txt",
                      "front_page_note": "/data/home/shannc/chula-stem/nextflow/reports/wes_front_page_note.txt"]


params.ref = ["genome": "/data/project/stemcell/shannc/reference/genomes/GRCh38.p14_filtered.fasta",
              "homopolymers_microsatellites": "/data/project/stemcell/shannc/reference/tool_specific/mishomopoly_GRCh38_filtered.tsv",
              "mappability": "/data/project/stemcell/shannc/reference/tool_specific/GCA_000001405.15_GRCh38_no_alt_analysis_set_100.bw",
              "targets": "/data/project/stemcell/shannc/reference/exome_kits/MGI_V5/MGI_Exome_Capture_V5_lifted_renamed.bed.gz",
              "baits": "/data/project/stemcell/shannc/reference/exome_kits/MGI_V5/MGI_Exome_Capture_V5_lifted_renamed.bed.gz", // Bait information not available
              "baits_il": "/data/project/stemcell/shannc/reference/exome_kits/SureSelectHumanAllExonV6Hg38/Covered.interval_list",
              "targets_il": "/data/project/stemcell/shannc/reference/exome_kits/MGI_V5/MGI_Exome_Capture_V5_lifted_renamed.interval_list"
              "baits_unzipped": "/data/project/stemcell/shannc/reference/exome_kits/MGI_V5/MGI_Exome_Capture_V5_lifted_renamed.bed.gz",
              "panel_of_normals": null,
              "snp_blacklist": null,
              "purecn_bait_intervals": null,
              "genome_blacklist": "/data/project/stemcell/shannc/reference/blacklists/ENCFF356LFX_renamed.bed",
              "delly_exclude": "",
              "known_variants": ["/data/project/stemcell/shannc/reference/variants/dbSNP_renamed_germline.vcf.gz", "/data/project/stemcell/shannc/reference/variants/gnomADv4.1.0_Exomes/random/gnomADv4.1_subset.vcf.gz"],
              "genome_gff": "/data/project/stemcell/shannc/reference/genomes/GRCh38.p14_filtered.gff",
              "add_to_pon": null,
              "germline": "/data/project/stemcell/shannc/reference/variants/gnomADv4.1.0_Exomes/biallelic/all.vcf.gz",
              "dbsnp": "/data/project/stemcell/shannc/reference/variants/dbSNP_renamed_germline.vcf.gz",
              "pileup": "/data/project/stemcell/shannc/reference/variants/gnomADv4.1.0_Exomes/biallelic/all.vcf.gz",
              "cnv_reference": "/data/project/stemcell/shannc/reference/variants/structural/aggregated_cnv.tsv",
              "msi_reference": "/data/project/stemcell/shannc/reference/variants/structural/aggregated_msi.tsv",
              "clingen_dosage": "/data/project/stemcell/shannc/reference/therapeutics/Clingen-Dosage-Sensitivity-2024-12-16.csv",
              "civic_cache": "/data/home/shannc/chula-stem/resources/civic_cache.json",
              "pandrugs2_cache": "/data/home/shannc/chula-stem/resources/2025-01-13_pandrugs2_cache_pdac.json",
              "cosmic_reference": "/data/home/shannc/chula-stem/resources/config/cosmic_signatures_v3.4-2024-12-26.csv",
              "clingen_gene": "/data/project/stemcell/shannc/reference/therapeutics/Clingen-Gene-Disease-Summary-2024-12-19.csv"]

// * Variant calling parameters

// ** Quality control filters
params.small_qc = ["accepted_filters": ["PASS"],
                   "min_tumor_depth": 0,
                   "max_normal_depth": 0,
                   "canonical": true,
                   "min_vaf": 0.10]
params.sv_qc = ["accepted_filters": ["PASS"]]
params.qc = ["accepted_filters": ["PASS"],
             "min_tumor_depth": 0,
             "max_normal_depth": 0,
             "min_vaf": 0.10,
             "min_callers": 2]

// * TEMPORARY
process {
    withName: "REPORT" {
        ext.misc = ["plot": false, "show_therapies": true]
        ext.extra_paths = ["therapy_data_cache": "/data/home/shannc/.cache/2025-01-13_therapy_data.parquet"]
    }
}

// * Misc
// Add gnomad subsets to known_variants
// gnomad_subsets = [
//     "$refdir/variants/gnomADv4.1.0_Exomes/random/chr10.vcf.bgz",
//     "$refdir/variants/gnomADv4.1.0_Exomes/random/chr11.vcf.bgz",
//     "$refdir/variants/gnomADv4.1.0_Exomes/random/chr12.vcf.bgz",
//     "$refdir/variants/gnomADv4.1.0_Exomes/random/chr13.vcf.bgz",
//     "$refdir/variants/gnomADv4.1.0_Exomes/random/chr14.vcf.bgz",
//     "$refdir/variants/gnomADv4.1.0_Exomes/random/chr15.vcf.bgz",
//     "$refdir/variants/gnomADv4.1.0_Exomes/random/chr16.vcf.bgz",
//     "$refdir/variants/gnomADv4.1.0_Exomes/random/chr17.vcf.bgz",
//     "$refdir/variants/gnomADv4.1.0_Exomes/random/chr18.vcf.bgz",
//     "$refdir/variants/gnomADv4.1.0_Exomes/random/chr19.vcf.bgz",
//     "$refdir/variants/gnomADv4.1.0_Exomes/random/chr1.vcf.bgz",
//     "$refdir/variants/gnomADv4.1.0_Exomes/random/chr20.vcf.bgz",
//     "$refdir/variants/gnomADv4.1.0_Exomes/random/chr21.vcf.bgz",
//     "$refdir/variants/gnomADv4.1.0_Exomes/random/chr22.vcf.bgz",
//     "$refdir/variants/gnomADv4.1.0_Exomes/random/chr2.vcf.bgz",
//     "$refdir/variants/gnomADv4.1.0_Exomes/random/chr3.vcf.bgz",
//     "$refdir/variants/gnomADv4.1.0_Exomes/random/chr4.vcf.bgz",
//     "$refdir/variants/gnomADv4.1.0_Exomes/random/chr5.vcf.bgz",
//     "$refdir/variants/gnomADv4.1.0_Exomes/random/chr6.vcf.bgz",
//     "$refdir/variants/gnomADv4.1.0_Exomes/random/chr7.vcf.bgz",
//     "$refdir/variants/gnomADv4.1.0_Exomes/random/chr8.vcf.bgz",
//     "$refdir/variants/gnomADv4.1.0_Exomes/random/chr9.vcf.bgz",
//     "$refdir/variants/gnomADv4.1.0_Exomes/random/chrX.vcf.bgz",
//     "$refdir/variants/gnomADv4.1.0_Exomes/random/chrY.vcf.bgz",
// ]
