process MARK_DUPLICATES {
    ext version: params.gatk_version
    

    publishDir "${meta.out}", mode:"copy", saveAs: params.saveFn
    publishDir "${meta.log}", mode:'copy', pattern: "*.log"

    input:
    tuple val(meta), path(aligned)
    val(reference)
    val(module_number)

    output:
    tuple val(meta), path(out), emit: dedup
    tuple val(meta), path(txt), emit: qc

    script:
    out = params.getName(module_number, meta, "dedup", "bam")
    txt = params.getName(module_number, meta, "dedup_metrics", "txt")
    check = file("${meta.out}/${out}")
    if (check.exists()) {
        """
        ln -sr ${check} .
        ln -sr ${meta.log}/dedup.log .
        ln -sr ${meta.out}/${txt} .
        """
    } else {
        """
        gatk SortSam -I $aligned \\
            -O sorted.bam \\
            --SORT_ORDER coordinate

        gatk ReorderSam -I sorted.bam \\
            -O reordered.bam \\
            --SEQUENCE_DICTIONARY $reference

        gatk MarkDuplicates \\
            -I sorted.bam \\
            --ASSUME_SORT_ORDER coordinate \\
            -M ${txt} \\
            -O ${out}

        get_nextflow_log.bash dedup.log
        """
    }

}
// GATK requires that reads are sorted in coordinate order
// Contigs are also expected to be ordered the same way across files
