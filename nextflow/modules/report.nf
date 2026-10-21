process REPORT {
    ext version: "1"
    maxForks 1 // Limit parallel execution to limit API calls

    input:
    tuple val(meta), val(path_map)
    val type

    output:

    //
    tuple val(meta), path(output), emit: report
    path "report_tmp"
    tuple val(meta), path("*.log"), emit: log

    script:
    //
    output = Utl.getName(meta, "Report", "pdf")
    paths = path_map.collectEntries { k, v ->
        [(k): v instanceof Path ? v.toString() : v]
    }
    all_meta = ["meta": meta + params.report_data, "paths": paths + task.ext.extra_paths, "misc": task.ext.misc, "other_text": params.report_text]
    specification = Utl.mapToJson(all_meta)
    """
    echo -e '${specification}' > spec.json

    pipeline_report -s spec.json \\
        -o ${output} \\
        -t ${type} \\
        -d ./report_tmp

    get_nextflow_log.bash report.log
    """
    //
}
