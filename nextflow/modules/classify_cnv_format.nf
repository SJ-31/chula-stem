process CLASSIFY_CNV_FORMAT {
    ext version: "1"

    
    

    input:
    tuple val(meta), path(cnv)
    
    //

    output:
    tuple val(meta), path(output), emit: bed
    tuple val(meta), path("*.log"), emit: log
    //

    script:
    caller = meta.caller
    output = Utl.getName(meta, "ClassifyCNV_format_${caller}", "bed")
    check = file("${meta.out}/${output}")
    args = task.ext.args.join(" ")
    if (check.exists()) {
        """
        ln -sr ${check} .
        ln -sr ${meta.log}/classify_cnv_format_${caller}.log .
        """
    } else {
        """
        classify_cnv_format -c ${caller} -i ${cnv} -o ${output}

        get_nextflow_log.bash classify_cnv_format_${caller}.log
        """
    }
    //
}
