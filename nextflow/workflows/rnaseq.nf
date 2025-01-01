include { PREPROCESS_FASTQ } from "../subworkflows/preprocess_fastq.nf"
include { KALLISTO } from '../modules/kallisto.nf'
include { STAR_FUSION } from "../modules/star_fusion.nf"

workflow RNASEQ {

    main:

    PREPROCESS_FASTQ(params.input, params.outdir, params.logdir, "rnaseq")

    // TODO: find a way to determine strandedness, per
    KALLISTO(PREPROCESS_FASTQ.out.trimmed, params.ref.kallisto_index, null, 2)
    STAR_FUSION(PREPROCESS_FASTQ.out.chimeric, params.ref.star_lib, 2)

}
