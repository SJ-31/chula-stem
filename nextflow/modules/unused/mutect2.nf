process MUTECT2 {
    ext version: params.gatk_version

    input:
    tuple val(meta), path(normal), path(tumor), path(indices, arity: "2")
    val reference
    val target_intervals // For exome data, path to target intervals file
    val germline_resource // A germline resource in the form of a population vcf with AF
    val panel_of_normals
    val interval_padding
    val tumor_only

    output:
    tuple val(meta), path(out), emit: variants
    tuple val(meta), path(raw), emit: raw
    tuple val(meta), path(stats), emit: stats
    tuple val(meta), path("*.log"), emit: logfile

    script:
    //
    out = Utl.getName(meta, "Mutect2", "vcf.gz")
    stats = "${out}.stats"
    raw = Utl.getName(meta, "Mutect2_raw", "tar.gz")
    target_flag = target_intervals ? " --intervals ${target_intervals} " : ""
    normal_flag = !tumor_only ? "-I ${normal} -normal ${meta.RGSM_normal} " : ""

    args = task.ext.args.join(" ")
    germline_flag = germline_resource ? " --germline-resource ${germline_resource} " : ""
    pon_flag = panel_of_normals ? " --panel-of-normals ${panel_of_normals}" : ""

    // Required for PureCN
    genotype_germline_flag = !tumor_only ? " --genotype-germline-sites true" : ""
    genotype_pon_flag = panel_of_normals ? "--genotype-pon-sites true" : ""
    """
    gatk Mutect2 \\
        ${args} \\
        -R ${reference} \\
        -I ${tumor} \\
        ${normal_flag} \\
        ${germline_flag} \\
        ${pon_flag} \\
        ${genotype_germline_flag} \\
        ${genotype_pon_flag} \\
        ${target_flag} \\
        --f1r2-tar-gz ${raw} \\
        --interval-padding ${interval_padding} \\
        --output tmp.vcf.gz

    mv tmp.vcf.gz.stats ${stats}

    if [[ "${tumor_only}" == "false" ]]; then
        bcftools view -s "${meta.RGSM_normal},${meta.RGSM_tumor}" tmp.vcf.gz | \\
            vcf_info_add_tag.bash -n ${params.source_name} \\
                -d "${params.source_description}" \\
                -b '.' \\
                -t String \\
                -a mutect2 \\
                -o ${out}
    else
        vcf_info_add_tag.bash -n ${params.source_name} \\
            -d "${params.source_description}" \\
            -b '.' \\
            -t String \\
            -a mutect2 \\
            -i tmp.vcf.gz \\
            -o ${out}
    fi

    get_nextflow_log.bash mutect2.log
    """
}
