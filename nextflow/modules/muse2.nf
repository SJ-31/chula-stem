process MUSE2 {
    ext version: "2.1.2"

    publishDir "$meta.out", mode:"copy", saveAs: params.saveFn
    publishDir "$meta.log", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(normal), path(tumor), path(indices, arity: "2")
    val(reference)
    val(dbsnp) // Docs ask for dbSNP vcf specifically
    val(omics_type)
    val(module_number)
    //

    output:
    tuple val(meta), path("${output}.gz"), emit: variants
    path("${prefix}.MuSE.txt")
    path("*.log")
    //

    script:
    prefix = "${module_number}-${meta.filename}"
    output = "${prefix}-MuSE.vcf"
    check = file("${meta.out}/${output}.gz")
    if (omics_type == "wgs") {
        data_flag = " -G "
    } else if (omics_type == "exome") {
        data_flag = " -E "
    } else {
        throw new Exception("Omics type must be either 'wgs' or 'exome'")
    }
    args = task.ext.args.join(" ")
    if (check.exists()) {
        """
        ln -sr ${check} .
        ln -sr ${meta.out}/${prefix}.MuSE.txt .
        ln -sr ${meta.log}/muse2.log .
        """
    } else {
        """
        MuSE call \\
            -f ${reference} \\
            -O ${prefix} \\
            -n ${task.ext.cores} \\
            ${tumor} \\
            ${normal}

        MuSE sump \\
            -I ${prefix}.MuSE.txt \\
            -O tmp.vcf \\
            -n ${task.ext.cores} \\
            -D ${dbsnp} \\
            ${data_flag}

        rename_vcf.bash -v -i tmp.vcf -o tmp2.vcf.gz \\
            -n "${meta.RGSM_normal}" -t "${meta.RGSM_tumor}"

        vcf_info_add_tag -n ${params.source_name} \\
            -d "$params.source_description" \\
            -b '.' \\
            -t String \\
            -a muse2 \\
            -i tmp2.vcf.gz \\
            -o ${output}

        bgzip ${output}
        get_nextflow_log.bash muse2.log
        """
    }
    //
}

