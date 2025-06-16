// * Required params
workDir = "/data/project/stemcell/shannc/output/nf-work/WES_TEST"
includeConfig "$projectDir/input/wes_tools.config.nf"
params.outdir = "/data/project/stemcell/shannc/output/WES_TEST"
params.input = "$projectDir/input/manifests/wes_test.csv"
// For debugging the pipeline ONLY, the actual
params.logdir = "${params.outdir}/log"
params.routine = "wes"
params.cohort = "test"
params.tumor_only = false

includeConfig "${projectDir}/input/references.config.nf"
params.report_data = ["hospital": "Dr.Nipan CU",
                      "physician": "Assoc. Prof. Nipan Israsena, M.D., Ph.D.",
                      "diagnosis": "test",
                      "sample_type": "organoid",
                      "collection_date": null,
]
params.report_text = ["front_page_details": "${projectDir}/reports/wes_front_page_details.txt",
                      "title": "Whole-exome sequencing results",
                      "disclaimer": "${projectDir}/reports/wes_disclaimer.txt",
                      "front_page_note": "${projectDir}/reports/wes_front_page_note.txt"]

params.ref = ["genome": genome,
              "homopolymers_microsatellites": homopolymers_microsatellites,
              "targets": targets,
              "baits": baits,
              "baits_il": baits_il,
              "targets_il": targets_il,
              "baits_unzipped": baits_unzipped,
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
params.qc = ["accepted_filters": ["PASS"],
             "min_tumor_depth": 10,
             "max_normal_depth": 10,
             "min_vaf": 0.10,
             "min_callers": 2,
             "informative": true,
             "canonical": true,
             "impact": true]

// * TEMPORARY
// <2025-01-08 Wed> Having difficulties with the therapy db caching, so for now
// will produce report without that data
process {
    withName: "REPORT" {
        ext.args = ["plot": false, "show_therapies": true]
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

