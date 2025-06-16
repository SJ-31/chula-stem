params.outdir = "/data/project/stemcell/shannc/tests/scRNASeq" // Directory for output
workDir = "$params.outdir/work"
params.input = "" // Path to manifest file of sample information. See routine-specific documentation for input format
params.logdir = "${params.outdir}/log" // Directory for logging process output
params.routine = "sc_rnaseq" // Name of routine to run
params.cohort = "test" // Optional: name of the sample cohort

refdir = "/data/project/stemcell/shannc/reference"
params.strandedness = ""
params.ref = [
    "star_index": "${refdir}/tool_specific/GRCh38_star_index",
    "barcodes": ,
    "genome_gff": "${refdir}/genomes/Homo_sapiens.GRCh38.113.gtf",
    "star_lib": "${refdir}/tool_specific/GRCh38_gencode_v44_CTAT_lib_Oct292023.plug-n-play.tar.gz"
    "genome": "${refdir}/genomes/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa",
    "transcriptome": "${refdir}/transcriptomes/Homo_sapiens.GRCh38.cdna.all.fa",
]

// * Variant calling parameters

// ** Quality control filters
params.small_qc = ["accepted_filters": ["PASS"],
                   "min_tumor_depth": 10,
                   "max_normal_depth": 10,
                   "min_VAF": 0.10]
params.sv_qc = ["accepted_filters": ["PASS"]]
