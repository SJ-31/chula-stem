process LEARN_READ_ORIENTATION {
    ext version: params.gatk_version

    input:
    tuple val(meta), path(raw)

    output:
    //
    tuple val(meta), path(output), emit: ro
    tuple val(meta), path("*.log"), emit: log

    script:
    //
    output = Utl.getName(meta, "RO_model", "tar.gz")
    """
    gatk LearnReadOrientationModel \\
        -I ${raw} \\
        -O ${output}

    get_nextflow_log.bash learn_read_orientation.log
    """
}
