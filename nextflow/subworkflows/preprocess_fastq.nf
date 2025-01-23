include { FASTP } from "../modules/fastp.nf"
include { BWA } from "../modules/bwa.nf"
include { MARK_DUPLICATES } from "../modules/mark_duplicates.nf"
include { BQSR } from "../modules/bqsr.nf"
include { SAMTOOLS_INDEX } from "../modules/samtools_index.nf"
include { STAR } from "../modules/star.nf"
include { STAR_SOLO } from "../modules/star_solo.nf"
include { SC_RNASEQ_QC } from "../modules/sc_rnaseq_qc.nf"

workflow PREPROCESS_FASTQ {
    take:
    manifest
    outdir
    logdir
    omics_type
    offset // Offset for module number

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

    FASTP(input, offset + 1)
    if (omics_type == "rnaseq") {
        STAR(FASTP.out.passed, params.ref.star_index, params.strandedness,
             "UniqueIdenticalNotMulti", true, params.ref.genome_gff, offset + 2)
        bams = STAR.out.mapped
        chimeric_junction = STAR.out.chimeric
        counts = STAR.out.counts
    } else if (omics_type == "sc_rnaseq") {
        STAR_SOLO(FASTP.out.passed, params.ref.star_index, params.strandedness, true,
                  params.ref.genome_gff, params.ref.barcodes, offset + 2)
        bams = STAR_SOLO.out.mapped
        chimeric_junction = STAR_SOLO.out.chimeric
        SC_RNASEQ_QC(STAR_SOLO.out.counts, params.ref.genome_gff_sqlite, offset + 3)
        counts = SC_RNASEQ_QC.out.counts
    } else {
        BWA(FASTP.out.passed, params.ref.genome, offset + 2)
        MARK_DUPLICATES(BWA.out.mapped, params.ref.genome, offset + 3)
        BQSR(MARK_DUPLICATES.out.dedup, params.ref.genome, params.ref.known_variants,
             false, offset + 4)
        bams = BQSR.out.bam
        chimeric_junction = Channel.empty()
        counts = Channel.empty()
    }
    SAMTOOLS_INDEX(bams) // indices are required by certain callers

    emit:
    bam = bams
    bam_index = SAMTOOLS_INDEX.out.index
    fastp_json = FASTP.out.json
    trimmed = FASTP.out.passed
    chimeric = chimeric_junction
    counts = counts
}
