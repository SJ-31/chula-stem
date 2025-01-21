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

    // Metrics
    //
    metric_bams = Utl.modifyMeta(PREPROCESS_FASTQ.out.bam,
                                 ["out": { "${params.outdir}/${it[0].id}/metrics" },
                                  "log": { "${params.logdir}/${it[0].id}/metrics" }])

    DUPRADAR(metric_bams, params.ref.genome_gff, params.strandedness, true, 4)

    to_picard = Utl.joinFirst(metric_bams, [PREPROCESS_FASTQ.out.bam_index])
    PICARD(to_picard, "rnaseq", params.ref.genome, "", "", params.ref.genome_gff, 4)

    to_multiqc = PREPROCESS_FASTQ.out.fastp_json.mix(PICARD.out.metrics,
                                                     PREPROCESS_FASTQ.out.counts
    )
        .flatten().collect().map({ [["out": params.outdir, "log": params.logdir,
                                     "filename": cohort_name], it] })
    MULTIQC(to_multiqc, 5)
}
