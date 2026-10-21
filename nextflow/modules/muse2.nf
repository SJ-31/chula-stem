process MUSE2 {
    ext version: "2.1.2"
    label "mid_mem"

    input:
    tuple val(meta), path(normal), path(tumor), path(indices, arity: "2")
    val reference
    val dbsnp // Docs ask for dbSNP vcf specifically
    val omics_type

    output:
    tuple val(meta), path(output), emit: variants
    tuple val(meta), path("${prefix}.MuSE.txt"), emit: txt
    tuple val(meta), path("*.log"), emit: log

    script:
    //
    prefix = Utl.getName(meta)
    output = "${prefix}-MuSE.vcf.gz"
    if (omics_type == "wgs") {
        data_flag = " -G "
    }
    else if (omics_type == "exome") {
        data_flag = " -E "
    }
    else {
        throw new Exception("Omics type must be either 'wgs' or 'exome'")
    }
    args = task.ext.args.join(" ")
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

    rename_vcf.bash -v -i tmp.vcf -n "${meta.RGSM_normal}" -t "${meta.RGSM_tumor}" | \\
        vcf_info_add_tag.bash -n ${params.source_name} \\
            -d "${params.source_description}" \\
            -b '.' \\
            -t String \\
            -a muse2 \\
            -o ${output}

    get_nextflow_log.bash muse2.log
    """
}
