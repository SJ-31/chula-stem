process DELLY_COV {
    ext version: "1.3.1"

    input:
    tuple val(meta), path(bam), path(indices, arity: "1..*")
    val reference
    val mappability

    output:

    //
    tuple val(meta), path(output), emit: cov
    tuple val(meta), path("*.log"), emit: logfile

    script:
    //
    output = Utl.getName(meta, "Delly", "cov.gz")

    """
    delly cnv -a \\
        --genome ${reference} \\
        --mappability ${mappability} \\
        --covfile ${output} \\
        --outfile tmp.bcf \\
        ${bam}

    get_nextflow_log.bash delly_cov.log
    """
}
