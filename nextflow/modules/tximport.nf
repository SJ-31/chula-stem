// Use the tximport package to aggregate transcript-level abundances in bulk
process TXIMPORT {
    ext version: "1"

    
    

    input:
    tuple val(meta), path(abundance)
    val(type)
    path(reference)
    
    //

    output:
    tuple val(meta), path(output)
    tuple val(meta), path("*.log"), emit: log
    //

    script:
    output = Utl.getName(meta, "imported", "tsv")
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
