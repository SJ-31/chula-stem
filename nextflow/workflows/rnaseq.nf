include { PREPROCESS_FASTQ } from "../subworkflows/preprocess_fastq.nf"
include { KALLISTO } from '../modules/kallisto.nf'
include { STAR_FUSION } from "../modules/star_fusion.nf"
include { MULTIQC } from "../modules/multiqc.nf"
include { DUPRADAR } from "../modules/dupradar.nf"
include { PICARD } from "../modules/picard.nf"

workflow rnaseq {

    main:

    def cohort_name = params.cohort ? params.cohort : "cohort"
    PREPROCESS_FASTQ(params.input, params.outdir, params.logdir, "rnaseq", 0)

    // params.strandedness is forward|reverse|unstranded
    // if (!params.strandedness) {
    //     // Insert TODO: https://github.com/signalbash/how_are_we_stranded_here
    // }
    KALLISTO(PREPROCESS_FASTQ.out.trimmed, params.ref.kallisto_index,
             params.strandedness, 2)
    if (params.detect_fusion) {
        STAR_FUSION(PREPROCESS_FASTQ.out.chimeric, params.ref.star_lib, 2)
    }
}
