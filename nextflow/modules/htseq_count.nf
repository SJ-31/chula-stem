process HTSEQ_COUNT {
    ext version: "2.0.5"

    input:
    tuple val(meta), path(alignments)
    val reference // Wants genome reference as gtf
    val sort_order // pos (position) or name
    val strandedness

    output:


    //
    tuple val(meta), path(output)
    tuple val(meta), path("*.log"), emit: logfile

    script:
    //
    output = Utl.getName(meta, "HTSeq_count", "tsv")
    if (strandedness == "forward") {
        strandedness_flag = "yes"
    }
    else if (strandedness == "reverse") {
        strandedness_flag = "reverse "
    }
    else {
        strandedness_flag = "no"
    }
    mode = task.ext.mode ? task.ext.mode : "union"
    args = task.ext.args.join(" ")
    """
    htseq-count ${alignments} ${reference} \\
        ${args} \\
        -r ${sort_order} \\
        -s ${strandedness_flag} \\
        -m ${mode} \\
        -c ${output} \\
        --with-header

    get_nextflow_log.bash htseq_count.log
    """
}
