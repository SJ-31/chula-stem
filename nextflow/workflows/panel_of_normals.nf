include { FASTP } from "../modules/fastp.nf"
include { BWA } from "../modules/bwa.nf"
include { MARK_DUPLICATES } from "../modules/mark_duplicates.nf"
include { BQSR } from "../modules/bqsr.nf"
include { SAMTOOLS_INDEX } from "../modules/samtools_index.nf"
include { MUTECT2 } from "../modules/mutect2.nf"
include { GENOMICS_DB_IMPORT } from '../modules/genomics_db_import.nf'
include { CREATE_PANEL_OF_NORMALS } from "../modules/create_panel_of_normals.nf"

process MUTECT2_TUMOR_ONLY {
    ext version: params.gatk_version

    publishDir "$meta.out", mode: "copy", saveAs: params.saveFn
    publishDir "$meta.log", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(normal), path(index)
    val(reference)
    val(target_intervals)
    val(germline)
    val(module_number)

    output:
    path(out)

    script:
    out = Utils.getName(module_number, meta, "Mutect2_normal", "vcf.gz")
    check = file("${meta.out}/${out}")
    target_flag = target_intervals != "" ? " --intervals ${target_intervals} " : ""
    if (check.exists()) {
        """
        ln -sr ${check} .
        """
    } else {
        """
        gatk Mutect2 -R ${reference} \\
            -I ${normal} \\
            -max-mnp-distance 0 \\
            ${target_flag} \\
            -O ${out}
        """
    }
}


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

    // TODO: modify empty to create a variable number of output
    // and also to take in paths, not vals

    empty_bams = EMPTY_FILES_1(BQSR.out.bam, 1).map(params.prependId)
    empty_indices = EMPTY_FILES_2(SAMTOOLS_INDEX.out.index).map(params.getId)

    to_mutect = empty_bams.join(BQSR.out.bam.map(params.getId))
        .join(SAMTOOLS_INDEX.out.index.map(params.getId))
        .join(empty_indices)
        .map({ it[1..-1] })

    // TODO: just use your existing mutect 2 for this
    MUTECT2(to_mutect, params.ref.genome, params.ref.targets,
                     params.ref.germline, 5)
    GENOMICS_DB_IMPORT(MUTECT2.out.somatic, params.ref.genome, params.ref.targets_il, 6)
    CREATE_PANEL_OF_NORMALS(GENOMICS_DB_IMPORT.out.db, params.ref.genome, 7)

}
