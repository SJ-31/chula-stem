// * Required params
workDir = "/data/project/stemcell/shannc/output/nf-work/HCC-RNASEQ"
includeConfig "$projectDir/input/rnaseq_tools.config.nf"
params.outdir = "/data/project/stemcell/shannc/output/HCC/RNASEQ"
if (params.debug) {
    params.input = "$projectDir/input/manifests/hcc_rnaseq_DEBUG.csv"
} else {
    params.input = "$projectDir/input/manifests/hcc_rnaseq_all.csv"
}
params.strandedness = "reverse" // From TruSeq's documentation, first read maps to antisense strand
params.relative_orientation = "inward"
params.detect_fusion = false
params.fastq_uncompress = "zcat" // Command used to uncompress fastq files for STAR
params.logdir = "${params.outdir}/log"
params.routine = "rnaseq"

refdir = "/data/project/stemcell/shannc/reference"

// <2025-01-09 Thu> Using Ensembl genome resources here for improved compatibility
// with Bioconductor packages
params.ref = [
    "genome": "${refdir}/genomes/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa",
    "transcriptome": "${refdir}/transcriptomes/Homo_sapiens.GRCh38.cdna.all.fa",
    "star_index": "${refdir}/tool_specific/GRCh38_star_index",
    "salmon_index": "${refdir}/tool_specific/salmon_index",
    "kallisto_index": "${refdir}/tool_specific/kallisto_GRCh38_ensembl.idx",
    "genome_gff": "${refdir}/genomes/Homo_sapiens.GRCh38.113.gtf",
    "genome_ref_flat": "${refdir}/genomes/Homo_sapiens.GRCh38.113.ref_flat",
    "star_lib": "${refdir}/tool_specific/GRCh38_gencode_v44_CTAT_lib_Oct292023.plug-n-play.tar.gz"
]
