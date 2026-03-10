process MUTECT2 {
    ext version: params.gatk_version

    publishDir "$meta.out", mode: "copy", saveAs: params.saveFn
    publishDir "$meta.log", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(normal), path(tumor), path(indices, arity: "2")
    val(reference)
    val(target_intervals) // For exome data, path to target intervals file
    val(germline_resource) // A germline resource in the form of a population vcf with AF
    val(panel_of_normals)
    val(interval_padding)
    val(module_number)
    //

    output:
    tuple val(meta), path(out), emit: variants
    tuple val(meta), path(raw), emit: raw
    tuple val(meta), path(stats), emit: stats
    path("*.log")
    //

    script:
    out = Utl.getName(module_number, meta, "Mutect2", "vcf.gz")
    stats = "${out}.stats"
    raw = Utl.getName(module_number, meta, "Mutect2_raw", "tar.gz")
    target_flag = target_intervals != "" ? " --intervals ${target_intervals} " : ""
    check = file("${meta.out}/${out}")
    normal_flag = !params.tumor_only ? "-I ${normal} -normal ${meta.RGSM_normal} " : ""

    args = task.ext.args.join(" ")
    germline_flag = germline_resource != "" ? " --germline-resource ${germline_resource} " : ""
    pon_flag = panel_of_normals ? " --panel-of-normals ${panel_of_normals}" : ""

    // Required for PureCN
    genotype_germline_flag = !params.tumor_only ? " --genotype-germline-sites true" : ""
    genotype_pon_flag = panel_of_normals ? "--genotype-pon-sites true" : ""

    if (check.exists()) {
        """
        ln -sr $check .
        ln -sr ${meta.log}/mutect2.log .
        ln -sr "${meta.out}/${stats}" .
        ln -sr "${meta.out}/${raw}" .
        """
    } else {
        """
        gatk Mutect2 \\
            ${args} \\
            -R $reference \\
            -I $tumor \\
            ${normal_flag} \\
            ${germline_flag} \\
            ${pon_flag} \\
            ${genotype_germline_flag} \\
            ${genotype_pon_flag} \\
            ${target_flag} \\
            --f1r2-tar-gz $raw \\
            --interval-padding ${interval_padding} \\
            --output tmp.vcf.gz

        mv tmp.vcf.gz.stats $stats

        if [[ "${params.tumor_only}" == "false" ]]; then
            bcftools view -s "${meta.RGSM_normal},${meta.RGSM_tumor}" tmp.vcf.gz | \\
                vcf_info_add_tag.bash -n ${params.source_name} \\
                    -d "$params.source_description" \\
                    -b '.' \\
                    -t String \\
                    -a mutect2 \\
                    -o ${out}
        else
            vcf_info_add_tag.bash -n ${params.source_name} \\
                -d "$params.source_description" \\
                -b '.' \\
                -t String \\
                -a mutect2 \\
                -i tmp.vcf.gz \\
                -o ${out}
        fi

        get_nextflow_log.bash mutect2.log
        """
    }
}
