include { PREPROCESS_FASTQ } from "../subworkflows/preprocess_fastq.nf"
include { KALLISTO } from '../modules/kallisto.nf'
include { STAR_FUSION } from "../modules/star_fusion.nf"

workflow rnaseq {

    main:

    PREPROCESS_FASTQ(params.input, params.outdir, params.logdir, "rnaseq")

    // params.strandedness is forward|reverse|unstranded
    // if (!params.strandedness) {
    //     // Insert TODO: https://github.com/signalbash/how_are_we_stranded_here
    // }
    KALLISTO(PREPROCESS_FASTQ.out.trimmed, params.ref.kallisto_index,
             params.strandedness, 2)
    STAR_FUSION(PREPROCESS_FASTQ.out.chimeric, params.ref.star_lib, 2)

}
