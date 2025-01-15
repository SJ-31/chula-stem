// Use the tximport package to aggregate transcript-level abundances in bulk
process TXIMPORT {
    ext version: "1"

    publishDir "${meta.out}", mode:"copy", saveAs: params.saveFn
    publishDir "${meta.log}", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(abundance)
    val(type)
    path(reference)
    val(module_number)
    //

    output:
    tuple val(meta), path(output)
    path("*.log")
    //

    script:
    output = Utl.getName(module_number, meta, "imported", "tsv")
    check = file("${meta.out}/${output}")
    args = task.ext.args.join(" ")
    if (check.exists()) {
        """
        ln -sr ${check} .
        ln -sr ${meta.log}/tximport.log .
        """
    } else {
        """
        tximport.R -i ${abundance} \\
            ${args} \\
            -o ${output} \\
            -r ${reference} \\
            -t ${type}

        get_nextflow_log.bash tximport.log
        """
    }
    //
}
