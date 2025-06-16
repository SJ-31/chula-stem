// * Required params
// Create a panel of normals to use for WES data
workDir = "/data/project/stemcell/shannc/output/nf-work/WES_PON"
includeConfig "$projectDir/input/wes_tumor_only_tools.config.nf"
params.outdir = "/data/project/stemcell/shannc/output/WES_PON"
params.input = "$projectDir/input/manifests/wes_normals.csv"
params.logdir = "${params.outdir}/log"
params.omics_type = "wes"
params.previous_vcfs = null
params.routine = "pon"
params.cohort = "all_normals"

process {
    withName: "MUTECT2" {
        ext.args = ["-max-mnp-distance 0"] // Recommended in GATK documentation
        // https://gatk.broadinstitute.org/hc/en-us/articles/360037227652-CreateSomaticPanelOfNormals-BETA
    }
}

includeConfig "${projectDir}/input/references.config.nf"

params.ref = ["genome": genome,
              "homopolymers_microsatellites": homopolymers_microsatellites,
              "targets": targets,
              "baits": baits,
              "baits_il": baits_il,
              "targets_il": targets_il,
              "baits_unzipped": baits_unzipped,
              "genome_blacklist": genome_blacklist,
              "known_variants": [dbsnp, gnomad_subset],
              "genome_gff": genome_gff,
              "germline": gnomad_subset_biallelic,
              "add_to_pon": [gnomad_subset_biallelic,
                             gatk_pon],
              "dbsnp": dbsnp,
              "pileup": gnomad_subset_biallelic,
              "cnv_reference": cnv_reference,
              "msi_reference": msi_reference,
              "clingen_dosage": clingen_dosage,
              "civic_cache": null,
              "pandrugs2_cache": null,
              "cosmic_reference": cosmic_reference,
              "clingen_gene": clingen_gene]

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
