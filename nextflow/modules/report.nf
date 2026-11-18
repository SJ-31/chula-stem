nextflow.enable.types = true

process REPORT {

    input:
    record(id: String,
           vep_small: Path,
           vcf_format_fields: Path,
           vep_sv: Path?,
           classify_cnv: Path?,
           concordant_cnvs: Path?,
           discordant_cnvs: Path?,
           mutational_signatures: Path?
        )
    report_type: String
    cache_path: String
    report_data: String

    output:

    //
    tuple val(meta), path(output), emit: report
    path "report_tmp"
    tuple val(meta), path("*.log"), emit: logfile

    script:
    //
    all_meta = ["meta":  + params.report_data, "paths": tool_outs, "misc": task.ext.misc, "other_text": params.report_text]
    specification = Utl.mapToJson(all_meta)
    """
    echo -e '' > spec.json

    pipeline_report -s spec.json \\
        -o ${output} \\
        -t ${type} \\
        -d ./report_tmp

    get_nextflow_log.bash report.log
    """
    //
}
