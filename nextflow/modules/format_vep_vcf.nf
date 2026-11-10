nextflow.enable.types = true

process FORMAT_VEP_VCF {
    input:
    record(
        id: String,
        vcf: Path,
    )
    fields: List<String>
    vep_field: String

    output:
    record(
        id: id,
        tsv: output,
    )

    script:
    field_arg = fields.collect{ v -> "-f ${v}" }.join(" ")
    output = "${vcf.simpleName}.tsv"
    """
    format_vep_vcf -i ${vcf} -o ${output} ${field_arg} -v ${vep_field}

    get_nextflow_log.bash format_vep_vcf.log
    """
}
