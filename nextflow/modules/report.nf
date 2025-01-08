process REPORT {
    ext version: "1"
    maxForks 1 // Limit parallel execution to limit API calls

    publishDir "${meta.out}", mode:"copy", saveAs: params.saveFn
    publishDir "${meta.log}", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), val(path_map)
    val(caches) // A list of caches to use
    val(type)
    val(module_number)
    //

    output:
    tuple val(meta), path(output)
    path("report_tmp")
    path("*.log")
    //

    script:
    output = Utl.getName(module_number, meta, "Report", "pdf")
    check = file("${meta.out}/${output}")
    paths = path_map.collectEntries( {
        k, v -> [(k): v instanceof Path ? v.toString() : v]
    } )
    all_meta = ["meta": meta, "paths": paths, "misc": task.ext.args,
                "other_text": params.report_text]
    specification = Utl.mapToJson(all_meta)
    """
    if [[ "${check.exists()}" == 'true' ]]; then
        cp -r ${meta.out}/report_tmp .
    fi

    echo -e '${specification}' > spec.json

    pipeline_report -s spec.json \\
        -o ${output} \\
        -t ${type} \\
        -d ./report_tmp

    get_nextflow_log.bash report.log
    """
    //
}
