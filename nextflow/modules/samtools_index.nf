process SAMTOOLS_INDEX {
    ext version: "1.21"

    input:
    tuple val(meta), path(bam)

    output:
    //
    tuple val(meta), path(output), emit: index
    tuple val(meta), path("*.log"), emit: log

    script:
    //
    output = "${bam}.bai"
    """
    samtools index ${bam}
    get_nextflow_log.bash samtools_index.log
    """
}
