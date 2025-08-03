// * Required params
params.tools = "" // Path to tool configuration
params.outdir = "" // Directory for output
params.input = "" // Path to manifest file of sample information. See routine-specific documentation for input format
params.logdir = "" // Directory for logging process output
params.routine = "" // Name of routine to run
params.cohort = "" // Optional: name of the sample cohort

// Map for pipeline resources
// See routine-specific documentation to see which need to be included
params.ref = []

// * Variant calling parameters

// ** Quality control filters
params.small_qc = ["accepted_filters": ["PASS"],
                   "min_tumor_depth": 10,
                   "max_normal_depth": 10,
                   "min_VAF": 0.10]
params.sv_qc = ["accepted_filters": ["PASS"]]
