process LEARN_READ_ORIENTATION {
    ext version: params.gatk_version

    publishDir "meta.out", mode:"copy", saveAs: { x -> x ==~ /.*\.log/ ? null : x }
    publishDir "meta.log", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(raw)
    val(module_number)
    //

    output:
    tuple val(meta), path(output), emit: ro
    path("*.log")
    //

    script:
    output = "${module_number}-${meta.filename}-RO_model.tar.gz"
    check = file("${meta.out}/${output}")
    args = task.ext.args.join(" ")
    if (check.exists()) {
        """
        ln -sr ${check} .
        ln -sr ${meta.log}/learn_read_orientation.log .
        """
    } else {
        """
        gatk LearnReadOrientationModel \\
            -I ${raw} \\
            -O ${output}

        get_nextflow_log.bash learn_read_orientation.log
        """
    }
    //
}
