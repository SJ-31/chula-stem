process MARK_DUPLICATES {
    publishDir "${meta.out}", mode:'copy'
    publishDir "${meta.log}", mode:'copy', pattern: "*.log"

    input:
    tuple val(meta), path(aligned)
    val(module_number)

    output:
    tuple val(meta), path("${meta.id}_dedup.bam"), emit: dedup
    tuple val(meta), path("${meta.id}_dedup_metrics.txt"), emit: qc

    script:
    out = "${module_number}-${meta.id}_dedup.bam"
    check = file("${meta.out}/${out}")
    if (check.exists()) {
        """
        cp $check.name .
        cp ${meta.log}/dedup.log .
        """
    } else {
        """
        samtools sort $aligned > sorted.bam
        gatk MarkDuplicatesSpark \\
            -I sorted.bam \\
            --ASSUME_SORT_ORDER coordinate \\
            -M ${meta.id}_dedup_metrics.txt \\
            -O ${out}

        cp .command.out dedup.log
        """
    }

}
