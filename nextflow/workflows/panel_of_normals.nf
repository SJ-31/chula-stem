include { PREPROCESS_FASTQ } from "../subworkflows/preprocess_fastq.nf"
include { MUTECT2 } from "../modules/mutect2.nf"
include { GENOMICS_DB_IMPORT } from '../modules/genomics_db_import.nf'
include { EMPTY_FILES as EMPTY_FILES_1 } from '../modules/empty_files.nf'
include { EMPTY_FILES as EMPTY_FILES_2 } from '../modules/empty_files.nf'
include { CLAIRS_TO } from "../modules/clairs_to.nf"
include { CREATE_PANEL_OF_NORMALS } from "../modules/create_panel_of_normals.nf"
include { OCTOPUS } from "../modules/octopus.nf"

params.tumor_only = true

workflow panel_of_normals {

    PREPROCESS_FASTQ(params.input, params.outdir, params.logdir, params.omics_type)
    if (params.previous_vcfs) {
        previous_vcfs = Channel.fromList(params.previous_vcfs.readLines()
                                        .collect({ file(it) }))
    } else {
        previous_vcfs = Channel.empty()
    }

    empty_bams = EMPTY_FILES_1(PREPROCESS_FASTQ.out.bam, 1)
    empty_indices = EMPTY_FILES_2(PREPROCESS_FASTQ.out.bam_index, 1)

    to_mutect2 = Utl.joinFirst(empty_bams, [PREPROCESS_FASTQ.out.bam,
                                         PREPROCESS_FASTQ.out.bam_index,
                                         empty_indices])
    to_clairs = Utl.joinFirst(PREPROCESS_FASTQ.out.bam,
                              [PREPROCESS_FASTQ.out.bam_index,
                               empty_indices])
    to_oct = Utl.joinFirst(to_mutect2, [empty_indices])

    MUTECT2(to_mutect2, params.ref.genome, params.ref.targets, params.ref.germline, 5)
    CLAIRS_TO(to_clairs, params.ref.genome, params.ref.targets, 5)
    OCTOPUS(to_oct, params.ref.genome, params.ref.targets, 5)

    all_variants = MUTECT2.out.variants.mix(CLAIRS_TO.out.variants, OCTOPUS.out.variants)
    to_genomics_db = all_variants.map({ it[1] }).mix(previous_vcfs).collect()
        .map({ [["filename": params.cohort,
                 "out": params.outdir,
                 "log": params.logdir ], it] })

    GENOMICS_DB_IMPORT(to_genomics_db, params.ref.genome, params.ref.targets_il, 6)
    CREATE_PANEL_OF_NORMALS(GENOMICS_DB_IMPORT.out.db, params.ref.genome, 7)

}
