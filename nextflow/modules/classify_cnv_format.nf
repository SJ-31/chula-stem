process CLASSIFY_CNV_FORMAT {
    ext version: "1"

    publishDir "${meta.out}", mode:"copy", saveAs: params.saveFn
    publishDir "${meta.log}", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(cnv)
    val(module_number)
    //

    output:
    tuple val(meta), path(output), emit: bed
    path("*.log")
    //

    script:
    caller = meta.caller
    output = Utl.getName(module_number, meta, "ClassifyCNV_format_${caller}", "bed")
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
