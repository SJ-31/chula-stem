include { FASTP } from "../modules/fastp.nf"
include { BWA } from "../modules/bwa.nf"
include { MARK_DUPLICATES } from "../modules/mark_duplicates.nf"
include { BQSR } from "../modules/bqsr.nf"
include { SAMTOOLS_INDEX } from "../modules/samtools_index.nf"
include { MUTECT2 } from "../modules/mutect2.nf"
include { GENOMICS_DB_IMPORT } from '../modules/genomics_db_import.nf'
include { EMPTY_FILES as EMPTY_FILES_1 } from '../modules/empty_files.nf'
include { EMPTY_FILES as EMPTY_FILES_2 } from '../modules/empty_files.nf'
include { CREATE_PANEL_OF_NORMALS } from "../modules/create_panel_of_normals.nf"

workflow panel_of_normals {

    main:
    input = Channel.fromPath(params.input)
        .splitCsv(header: true)
        .map { [["id": it.patient,
                "out": "${params.outdir}/${it.patient}/${it.source}",
                "type": it.source,
                "log": "${params.logdir}/${it.patient}/${it.source}",
                "filename": "${it.patient}_${it.source}", // Output filename
                "RGLB": it.read_group_library, // Read group info used by GATK tools
                "RGPL": it.read_group_platform,
                "RGPU": it.read_group_platform_unit,
                "RGSM": it.read_group_sample_name
            ], [file(it.fastq_1), file(it.fastq_2)]] }

    FASTP(input, 1)
    BWA(FASTP.out.passed, params.ref.genome, 2)
    MARK_DUPLICATES(BWA.out.mapped, params.ref.genome, 3)
    BQSR(MARK_DUPLICATES.out.dedup, params.ref.genome, params.ref.known_variants, 4)
    SAMTOOLS_INDEX(BQSR.out.bam)

    empty_bams = EMPTY_FILES_1(BQSR.out.bam, 1).map(params.prependId)
    empty_indices = EMPTY_FILES_2(SAMTOOLS_INDEX.out.index, 1).map(params.getId)

    to_mutect = empty_bams.join(BQSR.out.bam.map(params.getId))
        .join(SAMTOOLS_INDEX.out.index.map(params.getId))
        .join(empty_indices)
        .map({ it[1..-1] })

    MUTECT2(to_mutect, params.ref.genome, params.ref.targets, params.ref.germline, 5)
    GENOMICS_DB_IMPORT(MUTECT2.out.variants, params.ref.genome, params.ref.targets_il, 6)
    CREATE_PANEL_OF_NORMALS(GENOMICS_DB_IMPORT.out.db, params.ref.genome, 7)

}
