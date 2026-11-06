process MARK_DUPLICATES {
    ext version: params.gatk_version

    input:
    tuple val(meta), path(aligned)
    val reference

    output:
    tuple val(meta), path(out), emit: dedup
    tuple val(meta), path(txt), emit: qc

    script:
    out = Utl.getName(meta, "dedup", "bam")
    txt = Utl.getName(meta, "dedup_metrics", "txt")
    """
    gatk SortSam -I ${aligned} \\
        -O sorted.bam \\
        --SORT_ORDER coordinate

    gatk ReorderSam -I sorted.bam \\
        -O reordered.bam \\
        --SEQUENCE_DICTIONARY ${reference}

    gatk MarkDuplicates \\
        -I sorted.bam \\
        --ASSUME_SORT_ORDER coordinate \\
        -M ${txt} \\
        -O ${out}

    get_nextflow_log.bash dedup.log
    """
}
// GATK requires that reads are sorted in coordinate order
// Contigs are also expected to be ordered the same way across files
