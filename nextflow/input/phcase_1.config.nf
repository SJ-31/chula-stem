// * Required params
workDir = "/data/project/stemcell/shannc/output/nf-work/PHCase"
includeConfig "$projectDir/input/wes_tools.config.nf"
params.outdir = "/data/project/stemcell/shannc/output/PHCase"
params.input = "$projectDir/input/manifests/pdac_phcase1.csv"
params.logdir = "${params.outdir}/log"
params.routine = "wes"
params.cohort = "PHcase"
params.tumor_only = false
params.source_name = "SOURCE"

realdir = "/data/home/shannc/chula-stem/nextflow/"
includeConfig "${projectDir}/input/references.config.nf"
params.report_data = ["hospital": "Dr.Nipan CU",
                      "physician": "Assoc. Prof. Nipan Israsena, M.D., Ph.D.",
                      "diagnosis": "PDAC",
                      "sample_type": "organoid",
                      "collection_date": null,
]
params.report_text = ["front_page_details": "${realdir}/reports/wes_front_page_details.txt",
                      "title": "Whole-exome sequencing results",
                      "disclaimer": "${realdir}/reports/wes_disclaimer.txt",
                      "front_page_note": "${realdir}/reports/wes_front_page_note.txt"]

mgi_kit = "$refdir/exome_kits/MGI_V5"
mgi_targets_unzipped = "$mgi_kit/MGI_Exome_Capture_V5_lifted_renamed.bed"
mgi_targets = "$mgi_kit/MGI_Exome_Capture_V5_lifted_renamed.bed.gz"
mgi_targets_il = "$mgi_kit/MGI_Exome_Capture_V5_lifted_renamed.interval_list"
resourcedir = "/data/home/shannc/chula-stem/resources"

params.ref = ["genome": genome,
              "homopolymers_microsatellites": homopolymers_microsatellites,
              "mappability": "${refdir}/tool_specific/GCA_000001405.15_GRCh38_no_alt_analysis_set_100.bw",
              "targets": mgi_targets,
              "baits": mgi_targets, // Bait information not available
              "baits_il": mgi_targets_il,
              "targets_il": mgi_targets_il,
              "baits_unzipped": mgi_targets_unzipped,
              "panel_of_normals": null,
              "snp_blacklist": null,
              "purecn_bait_intervals": null,
              "genome_blacklist": genome_blacklist,
              "delly_exclude": "",
              "known_variants": [dbsnp, gnomad_subset],
              "genome_gff": genome_gff,
              "germline": gnomad_subset_biallelic,
              "dbsnp": dbsnp,
              "pileup": gnomad_subset_biallelic,
              "cnv_reference": cnv_reference,
              "msi_reference": msi_reference,
              "clingen_dosage": clingen_dosage,
              "civic_cache": "${resourcedir}/civic_cache.json",
              "pandrugs2_cache": "${resourcedir}/2025-01-13_pandrugs2_cache_pdac.json",
              "cosmic_reference": cosmic_reference,
              "clingen_gene": clingen_gene]

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

