nextflow.enable.types = true

process FORMAT_VEP_VCF {
    input:
    record(
        id: String,
        vcf: Path,
    )

    output:
    record(
        id: id,
        tsv: output,
    )

    script:
    output = "${vcf.simpleName}.tsv"
    """
    format_vep_vcf -i ${vcf} -o ${output} -f CALLERS -f NCALLERS -f SVTYPE

    get_nextflow_log.bash format_vep_vcf.log
    """
}
