nextflow.enable.types = true

process EXTRACT_FORMAT_FIELDS {
    input:
    record(
        id: String,
        normal_sample: String,
        tumor_sample: String,
        caller2vcfs: Map<String, Path>,
        suffix: String
    )

    output:
    record(
        id: id,
        tsv: output,
    )

    topic:
    record(id: id, process: "extract_format_fields", file: file("extract_format_fields.log")) >> "logs"

    script:
    //
    input_str = caller2vcfs.entrySet().collect { e -> "-i ${e.key}:${e.value}" }.join(" ")
    suffix = suffix ?: ""
    samples_field = ""
    if (normal_sample) {
        samples_field = "${samples_field} -s NORMAL -s ${id}_${normal_sample} -r NORMAL=${normal_sample}"
    }
    if (tumor_sample) {
        samples_field = "${samples_field} -s TUMOR -s ${id}_${tumor_sample} -r TUMOR=${id}_${tumor_sample}"
    }
    output = "${id}-format_fields${suffix}.tsv"
    """
    extract_vcf_format ${input_str} -o ${output} ${samples_field} 

    get_nextflow_log.bash extract_format_fields.log
    """
}
