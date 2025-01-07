process REPORT {
    ext version: "1"
    maxForks 1 // Limit parallel execution to limit API calls

    publishDir "${meta.out}", mode:"copy", saveAs: params.saveFn
    publishDir "${meta.log}", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta_map), val(path_map)
    val(type)
    val(module_number)
    //

    output:
    tuple val(meta_map), path(output)
    path("report_tmp")
    path("*.log")
    //

    script:
    output = Utl.getName(module_number, meta_map, "Report", "pdf")
    check = file("${meta_map.out}/${output}")
    if (check.exists()) {
        """
        ln -sr ${check} .
        cp -r ${meta_map.out}/report_tmp .
        ln -sr ${meta_map.log}/report.log .
        """
    } else {
        paths = path_map.collectEntries( {
            k, v -> [(k): v instanceof Path ? v.toString() : v]
        } )
        all_meta = [meta: meta_map,
                    paths: paths, misc: task.ext.args]
        specification = Utl.mapToJson(all_meta)
        """
        echo -e '${specification}' > spec.json

        pipeline_report -s spec.json \\
            -o ${output} \\
            -t ${type} \\
            -d ./report_tmp

        get_nextflow_log.bash report.log
        """
    }
    //
}
