process GRIDSS {
    ext version: "2.13.2"
    errorStrategy "ignore" // <2025-03-05 Wed> bug: gridss wasn't formatting
    // the vcf file correctly

    label "big_mem"

    input:
    tuple val(meta), path(normal), path(tumor), path(indices, arity: "2")
    val reference
    val blacklist // list of genomic regions to exclude, gridss docs recommends using

    output:
    // ENCODE blacklists

    //
    tuple val(meta), path(filtered), emit: variants
    path all
    tuple val(meta), path("*.log"), emit: logfile

    script:
    //
    all = Utl.getName(meta, "Gridss_all", "vcf.gz")
    filtered = Utl.getName(meta, "Gridss_confident", "vcf.gz")
    bam_input = !params.tumor_only ? " ${normal} ${tumor} " : tumor
    """
    gridss \\
        --threads ${task.cpus} \\
        --reference ${reference} \\
        --output tmp.vcf \\
        ${bam_input}

    gridss_somatic_filter \\
        --input tmp.vcf \\
        --output confident.vcf \\
        --fulloutput all.vcf

    names=("${filtered}" "${all}")
    files=(confident.vcf.bgz all.vcf.bgz)

    for i in \$(seq 0 1); do
        name="\${names[i]}"
        file="\${files[i]}"

        if [[ "\${file}" == "confident.vcf.bgz" && ${params.tumor_only} == "false" ]]; then
            bcftools view -s "${meta.RGSM_normal},${meta.RGSM_tumor}" -O z "\${file}" > tmp.vcf.gz
            f="tmp.vcf.gz"
        else
            f="\${file}"
        fi

        vcf_info_add_tag.bash -n ${params.source_name} \\
            -d '${params.source_description}' \\
            -b '.' \\
            -t String \\
            -a gridss \\
            -i "\${f}" \\
            -o "\$name"
    done

    get_nextflow_log.bash gridss.log
    """

    // TOOD: you can supply a pon with this
}
