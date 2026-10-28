include { PREPROCESS_FASTQ } from "../subworkflows/preprocess_fastq.nf"

// TODO: unfinished
//
workflow sc_rnaseq {
    // preprocess and align (to whole genome based on recommendations)
    // Can probably use fastp with less stringent flags
    // platform is 10x genomics, so will need to deal with cell barcodes somehow
    // Ask about UMI sequences

    // STAR_SOLO carries out cell counting and umi resolution,
    //

    //
    // IMPORTANT: modify ext.args for STAR solo to ensure that the CBs and UMIs
    // are all correct
    PREPROCESS_FASTQ(params.input, params.outdir, params.logdir, "sc_rnaseq", 0)

}
