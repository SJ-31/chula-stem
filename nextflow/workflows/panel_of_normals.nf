include { PREPROCESS_FASTQ } from "../subworkflows/preprocess_fastq.nf"
include { PANEL_OF_NORMALS_FROM_BAM } from "../subworkflows/panel_of_normals_from_bam.nf"

workflow panel_of_normals {

    PREPROCESS_FASTQ(params.input, params.outdir, params.logdir, params.omics_type, 0)
    if (params.previous_vcfs) {
        previous_vcfs = channel.fromList(params.previous_vcfs.readLines()
                                        .collect({ it -> file(it) }))
    } else {
        previous_vcfs = channel.empty()
    }
    def cohort_name = params.cohort ? params.cohort : "cohort"

    PANEL_OF_NORMALS_FROM_BAM(
        PREPROCESS_FASTQ.out.bam,
        PREPROCESS_FASTQ.out.bam_index,
        cohort_name,
        previous_vcfs,
        true
    )
}
