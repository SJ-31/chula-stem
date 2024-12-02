include { MUTECT2 } from "../modules/mutect2.nf"
include { LEARN_READ_ORIENTATION } from "../modules/learn_read_orientation.nf"
include { GET_PILEUP_SUMMARIES } from "../modules/get_pileup_summaries.nf"
include { CALCULATE_CONTAMINATION } from "../modules/calculate_contamination.nf"
include { FILTER_MUTECT_CALLS } from "../modules/filter_mutect_calls.nf"

// Subworkflow for somatic variant calling with Gatk's mutect2 and filtering its output
// according to best practices
workflow MUTECT2_COMPLETE {
    take:
    meta_and_bam
    module_number

    main:
    MUTECT2(meta_and_bam, params.ref.genome, params.ref.targets, module_number)
    LEARN_READ_ORIENTATION(MUTECT2.out.raw, module_number)

    to_pileup = meta_and_bam.map({ [it[0], it[2], it[3] ] })

    GET_PILEUP_SUMMARIES(to_pileup, params.ref.pileup, module_number)

    CALCULATE_CONTAMINATION(GET_PILEUP_SUMMARIES.out.pileup, module_number)

    to_filter = MUTECT2.out.variants.map({ [it[0].id] + it })
        .join(LEARN_READ_ORIENTATION.out.ro.map(params.getId))
        .join(CALCULATE_CONTAMINATION.out.c.map(params.getId))
        .join(CALCULATE_CONTAMINATION.out.s.map(params.getId))
        .join(MUTECT2.out.stats.map(params.getId))
        .map { it[1..-1] }

    FILTER_MUTECT_CALLS(to_filter, params.ref.genome, module_number)

    emit:
    FILTER_MUTECT_CALLS.out.filtered

}
