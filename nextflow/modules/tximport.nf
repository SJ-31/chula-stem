// Use the tximport package to aggregate transcript-level abundances in bulk
process TXIMPORT {
    ext version: "1"

    input:
    tuple val(meta), path(abundance)
    val type
    path reference

    output:

    //
    tuple val(meta), path(output)
    tuple val(meta), path("*.log"), emit: logfile

    script:
    //
    output = Utl.getName(meta, "imported", "tsv")
    args = task.ext.args.join(" ")

    """
    tximport.R -i ${abundance} \\
        ${args} \\
        -o ${output} \\
        -r ${reference} \\
        -t ${type}

    get_nextflow_log.bash tximport.log
    """
}
