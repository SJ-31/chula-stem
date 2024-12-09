include { FASTP } from "../modules/fastp.nf"
include { BWA } from "../modules/bwa.nf"
include { MARK_DUPLICATES } from "../modules/mark_duplicates.nf"
include { BQSR } from "../modules/bqsr.nf"
include { SAMTOOLS_INDEX } from "../modules/samtools_index.nf"

workflow  PREPROCESS_FASTQ {
    take:
    manifest
    outdir
    logdir

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
    BWA(FASTP.out.passed, params.ref.genome, 2)
    MARK_DUPLICATES(BWA.out.mapped, params.ref.genome, 3)
    BQSR(MARK_DUPLICATES.out.dedup, params.ref.genome, params.ref.known_variants, 4)
    SAMTOOLS_INDEX(BQSR.out.bam) // indices are required by certain callers


    emit:
    bam = BQSR.out.bam
    bam_index = SAMTOOLS_INDEX.out.index
    fastp_json = FASTP.out.json

}
