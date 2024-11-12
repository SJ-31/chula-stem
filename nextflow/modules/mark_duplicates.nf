process MARK_DUPLICATES {
    ext version: "4.6.1.0"
    

    publishDir "${meta.out}", mode:'copy'
    publishDir "${meta.log}", mode:'copy', pattern: "*.log"

    input:
    tuple val(meta), path(aligned)
    val(reference)
    val(module_number)

    output:
    tuple val(meta), path("${meta.id}_dedup.bam"), emit: dedup
    tuple val(meta), path("${meta.id}_dedup_metrics.txt"), emit: qc

    script:
    out = "${module_number}-${meta.id}_dedup.bam"
    check = file("${meta.out}/${out}")
    if (check.exists()) {
        """
        cp $check .
        cp ${meta.log}/dedup.log .
        """
    } else {
        """
        gatk SortSam -I $aligned \\
            -O sorted.bam \\
            --SORT_ORDER coordinate

        gatk ReorderSam -I sorted.bam \\
            -O reordered.bam \\
            --SEQUENCE_DICTIONARY $reference

        gatk MarkDuplicatesSpark \\
            -I sorted.bam \\
            --ASSUME_SORT_ORDER coordinate \\
            -M ${meta.id}_dedup_metrics.txt \\
            -O ${out}

        cp .command.out dedup.log
        """
    }

}
// GATK requires that reads are sorted in coordinate order
// Contigs are also expected to be ordered the same way across files
