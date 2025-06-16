// * Required params
includeConfig "$projectDir/input/pon_tools.groovy" // Path to tool configuration
params.outdir = "/data/project/stemcell/shannc/output/panel_of_normals" // Directory for output
params.input = "" // Path to manifest file of sample information. See routine-specific documentation for input format
params.logdir = "/data/project/stemcell/shannc/output/panel_of_normals/log" // Directory for logging process output
params.routine = "pon" // Name of routine to run
params.cohort = "cohort" // Optional: name of the sample cohort
includeConfig includeConfig // configuration for tools specific to the current routine

// Map for pipeline resources
// See routine-specific documentation to see which need to be included
refdir = "/data/project/stemcell/shannc/reference"
params.ref = ["genome": "$refdir/genomes/GRCh38.p14_filtered.fna",
              "baits_il": "$refdir/exome_kits/SureSelectHumanAllExonV6Hg38/Covered.interval_list",
              "targets_il": "$refdir/exome_kits/SureSelectHumanAllExonV6Hg38/Regions.interval_list",
              "baits_unzipped": "$refdir/exome_kits/SureSelectHumanAllExonV6Hg38/Unzipped_covered.bed",
              "known_variants": ["$refdir/variants/dbSNP_renamed_germline.vcf.gz",
                                 "$refdir/variants/gnomADv4.1.0_all.vcf.gz"],
              "dbsnp": "$refdir/variants/dbSNP_renamed_germline.vcf.gz",
]


// * Variant calling parameters

// ** Quality control filters
params.small_qc = ["accepted_filters": ["PASS"],
                   "min_tumor_depth": 10,
                   "max_normal_depth": 10,
                   "min_VAF": 0.10]
params.sv_qc = ["accepted_filters": ["PASS"]]
