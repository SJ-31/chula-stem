include { PREPROCESS_FASTQ } from "../subworkflows/preprocess_fastq.nf"
include { MUTECT2 } from "../modules/mutect2.nf"
include { GENOMICS_DB_IMPORT } from '../modules/genomics_db_import.nf'
include { EMPTY_FILES as EMPTY_FILES_1 } from '../modules/empty_files.nf'
include { EMPTY_FILES as EMPTY_FILES_2 } from '../modules/empty_files.nf'
include { EMPTY_FILES as EMPTY_FILES_3 } from '../modules/empty_files.nf'
include { CLAIRS_TO } from "../modules/clairs_to.nf"
include { CREATE_PANEL_OF_NORMALS } from "../modules/create_panel_of_normals.nf"
include { OCTOPUS } from "../modules/octopus.nf"
include { CNVKIT_PREP } from "../modules/cnvkit_prep.nf"
include { CONCAT_VCF } from "../modules/concat_vcf.nf"

workflow panel_of_normals {

    PREPROCESS_FASTQ(params.input, params.outdir, params.logdir, params.omics_type, 0)
    if (params.previous_vcfs) {
        previous_vcfs = Channel.fromList(params.previous_vcfs.readLines()
                                        .collect({ file(it) }))
    } else {
        previous_vcfs = Channel.empty()
    }
    def cohort_name = params.cohort ? params.cohort : "cohort"

    empty_bams = EMPTY_FILES_1(PREPROCESS_FASTQ.out.bam, "", 1)
    empty_indices_1 = EMPTY_FILES_2(PREPROCESS_FASTQ.out.bam_index, "E1", 1)
    empty_indices_2 = EMPTY_FILES_3(PREPROCESS_FASTQ.out.bam_index, "E2", 1)

    collected_normals = PREPROCESS_FASTQ.out.bam_index.mix(PREPROCESS_FASTQ.out.bam)
        .map({ it[1] }).toList()
    CNVKIT_PREP(Channel.of(["filename": cohort_name,
                            "out": "${params.outdir}/cnvkit_reference",
                            "log": "${params.outdir}/cnvkit_reference"])
                .merge(collected_normals) { meta, bams -> tuple(meta, tuple(bams)) },
                params.ref.genome, params.ref.baits_unzipped,
                params.ref.genome_blacklist, true, "hybrid", 4)

    to_mutect2 = Utl.joinFirst(empty_bams, [PREPROCESS_FASTQ.out.bam,
                                            PREPROCESS_FASTQ.out.bam_index,
                                            empty_indices_1])
        .map({it[0..-3] + [it[-2..-1]]})
    to_clairs = Utl.joinFirst(PREPROCESS_FASTQ.out.bam,
                              [PREPROCESS_FASTQ.out.bam_index,
                               empty_indices_1])
    to_oct = Utl.joinFirst(to_mutect2, [empty_indices_2])

    MUTECT2(to_mutect2, params.ref.genome, params.ref.targets, params.ref.germline, 5)
    CLAIRS_TO(to_clairs, params.ref.genome, params.ref.targets, 5)
    OCTOPUS(to_oct, params.ref.genome, params.ref.targets, 5)

    to_concat = Utl.joinFirst(MUTECT2.out.variants,
                              [CLAIRS_TO.out.variants, OCTOPUS.out.variants])
        .map({ [it[0], it[1..-1]] })

    CONCAT_VCF(to_concat, 5)

    // BUG: <2025-01-22 Wed> db_import isn't working, will try see if it's due to
    // handling the vcfs of the other callers
    vcf_channels = CONCAT_VCF.out.vcf.map({ it[1] }).mix(previous_vcfs).collect()

    // vcf_channels = MUTECT2.out.variants.map({ it[1] }).collect()
    to_pon = vcf_channels.map({ [["filename": params.cohort,
                                          "out": params.outdir,
                                          "log": params.logdir ], it] })

    min_samples = params.minimum_samples ? params.minimum_samples : 2
    CREATE_PANEL_OF_NORMALS(to_pon, min_samples, params.ref.add_to_pon, 7)
}
