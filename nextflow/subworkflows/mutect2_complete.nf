include { MUTECT2 } from "../modules/mutect2.nf"
include { LEARN_READ_ORIENTATION } from "../modules/learn_read_orientation.nf"
include { GET_PILEUP_SUMMARIES } from "../modules/get_pileup_summaries.nf"
include { CALCULATE_CONTAMINATION } from "../modules/calculate_contamination.nf"
include { FILTER_MUTECT_CALLS } from "../modules/filter_mutect_calls.nf"

// Subworkflow for somatic variant calling with Gatk's mutect2 and filtering its output
// according to best practices
// See https://gatk.broadinstitute.org/hc/en-us/articles/360035531132--How-to-Call-somatic-mutations-using-GATK4-Mutect2
workflow MUTECT2_COMPLETE {
    take:
    meta_and_bam
    panel_of_normals
    module_number

    main:
    MUTECT2(meta_and_bam,
            params.ref.genome,
            params.ref.targets,
            params.ref.germline,
            panel_of_normals,
            params.interval_padding ? params.interval_padding : 0,
            params.tumor_only,
            module_number)
    LEARN_READ_ORIENTATION(MUTECT2.out.raw, module_number)

    to_pileup = meta_and_bam.map({ it -> [it[0], it[2], it[3] ] })

    GET_PILEUP_SUMMARIES(to_pileup, params.ref.pileup, module_number)

    CALCULATE_CONTAMINATION(GET_PILEUP_SUMMARIES.out.pileup, module_number)

    to_filter = Utl.joinFirst(MUTECT2.out.variants,
                              [LEARN_READ_ORIENTATION.out.ro,
                               CALCULATE_CONTAMINATION.out.c,
                               CALCULATE_CONTAMINATION.out.s,
                               MUTECT2.out.stats])

    FILTER_MUTECT_CALLS(to_filter, params.ref.genome, module_number)

    emit:
    unfiltered = MUTECT2.out.variants
    filtered = FILTER_MUTECT_CALLS.out.filtered

}
