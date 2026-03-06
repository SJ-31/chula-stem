include { PREPROCESS_FASTQ } from "../subworkflows/preprocess_fastq.nf"
include { KALLISTO } from '../modules/kallisto.nf'
include { STAR_FUSION } from "../modules/star_fusion.nf"
include { MULTIQC } from "../modules/multiqc.nf"
include { DUPRADAR } from "../modules/dupradar.nf"
include { PICARD } from "../modules/picard.nf"
include { HTSEQ_COUNT } from "../modules/htseq_count.nf"
include { SALMON } from "../modules/salmon.nf"
include { COMBINE_COUNTS } from "../modules/combine_counts.nf"

workflow rnaseq {

    main:

    def cohort_name = params.cohort ? params.cohort : "cohort"
    PREPROCESS_FASTQ(params.input, params.outdir, params.logdir, "rnaseq", 0)

    // params.strandedness is forward|reverse|unstranded
    KALLISTO(PREPROCESS_FASTQ.out.trimmed, params.ref.kallisto_index,
             params.strandedness, 2)
    HTSEQ_COUNT(PREPROCESS_FASTQ.out.bam, params.ref.genome_gff, "pos",
                params.strandedness, 2)

    // params.relative_orientation is inward|outward|matching
    // SALMON(PREPROCESS_FASTQ.out.trimmed, params.ref.salmon_index, "",
    //        params.genome_gff, params.strandedness, params.relative_orientation, 2)
    //        BUG: salmon does not recognize files

    if (params.detect_fusion) {
        STAR_FUSION(PREPROCESS_FASTQ.out.chimeric, params.ref.star_lib, 2)
    }


    to_combine = PREPROCESS_FASTQ.out.counts.map({ it -> ["sample": it[0].id,
                                                    "type": it[0].type,
                                                    "file": it[1]] })
        .collect().map({ it -> [["out": params.outdir, "log": params.logdir,
                                 "filename": cohort_name], it] })
    COMBINE_COUNTS(to_combine, 4)
    
    

    // Metrics
    //
    metric_bams = Utl.modifyMeta(PREPROCESS_FASTQ.out.bam,
                                 ["out": { it -> "${params.outdir}/${it.id}/metrics" },
                                  "log": { it -> "${params.logdir}/${it.id}/metrics" }])

    DUPRADAR(metric_bams, params.ref.genome_gff, params.strandedness, true, 4)

    to_picard = Utl.joinFirst(metric_bams, [PREPROCESS_FASTQ.out.bam_index])
    PICARD(to_picard, "rnaseq", params.ref.genome, "", "", params.strandedness,
           params.ref.genome_ref_flat, 4)

    
    to_multiqc = PREPROCESS_FASTQ.out.fastp_json.mix(PICARD.out.metrics)
        .flatten().collect().map({ it -> [["out": params.outdir,
                                     "log": params.logdir,
                                     "filename": cohort_name], it] })
    MULTIQC(to_multiqc, 5)
}
