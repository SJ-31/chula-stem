include { FASTP } from "../modules/fastp.nf"
include { BWA } from "../modules/bwa.nf"
include { MARK_DUPLICATES } from "../modules/mark_duplicates.nf"
include { BQSR } from "../modules/bqsr.nf"
include { SAMTOOLS_INDEX } from "../modules/samtools_index.nf"
include { STAR } from "../modules/star.nf"

workflow  PREPROCESS_FASTQ {
    take:
    manifest
    outdir
    logdir
    data_type

    main:
    input = Channel.fromPath(manifest)
        .splitCsv(header: true)
        .map { [["id": it.patient,
                "out": "${outdir}/${it.patient}/${it.source}",
                "type": it.source,
                "log": "${logdir}/${it.patient}/${it.source}",
                "filename": "${it.patient}_${it.source}", // Output filename
                "RGLB": it.read_group_library, // Read group info used by GATK tools
                "RGPL": it.read_group_platform,
                "RGPU": it.read_group_platform_unit,
                "RGSM": it.read_group_sample_name
            ], [file(it.fastq_1), file(it.fastq_2)]] }

    FASTP(input, 1)
    if (data_type == "rnaseq") {
        STAR(FASTP.out.passed, params.ref.star_index, params.strandedness, 2)
        MARK_DUPLICATES(STAR.out.mapped, params.ref.genome, 3)
        bams = MARK_DUPLICATES.out.dedup
        chimeric_junction = STAR.out.chimeric
    } else {
        BWA(FASTP.out.passed, params.ref.genome, 2)
        MARK_DUPLICATES(BWA.out.mapped, params.ref.genome, 3)
        BQSR(MARK_DUPLICATES.out.dedup, params.ref.genome, params.ref.known_variants,
             false, 4)
        bams = BQSR.out.bam
        chimeric_junction = Channel.empty()
    }
    SAMTOOLS_INDEX(bams) // indices are required by certain callers

    emit:
    bam = bams
    bam_index = SAMTOOLS_INDEX.out.index
    fastp_json = FASTP.out.json
    trimmed = FASTP.out.passed
    chimeric = chimeric_junction
}
