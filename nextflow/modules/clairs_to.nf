process CLAIRS_TO {
    ext version: "0.3.1"

    label "big_mem"

    input:
    tuple val(meta), path(tumor), path(index), path(previous_variants)
    val reference
    val target_intervals

    output:

    //
    tuple val(meta), path(output), emit: variants
    tuple val(meta), path("*.log"), emit: log

    script:
    //
    output = Utl.getName(meta, "ClairS-TO", "vcf.gz")
    target_flag = target_intervals != "" ? " --bed_fn ${target_intervals} " : ""
    prev_flag = previous_variants.size() != 0 ? "--hybrid_mode_vcf_fn ${previous_variants}" : ""
    args = task.ext.args.join(" ")
    """
    name=\$( tag_from_sam.bash ${tumor} RG SM )
    run_clairs_to ${args} \\
        --tumor_bam_fn ${tumor} \\
        --ref_fn ${reference} \\
        --sample_name \${name} \\
        --threads ${task.cpus} \\
        --platform ${task.ext.platform} \\
        ${target_flag} \\
        ${prev_flag} \\
        --output_dir .

    if [[ -e snv.vcf.gz && -e indel.vcf.gz ]]; then
        bcftools concat snv.vcf.gz indel.vcf.gz -a | \\
            vcf_info_add_tag.bash -n ${params.source_name} \\
                -d "${params.source_description}" \\
                -b '.' \\
                -t String \\
                -a clairs-to \\
                -o tmp.vcf.gz
    else
        if [[ -e snv.vcf.gz ]]; then
            result=snv.vcf.gz
        else
            result=indel.vcf.gz
        fi
        vcf_info_add_tag.bash -n ${params.source_name} \\
            -d "${params.source_description}" \\
            -b '.' \\
            -i \$result \\
            -t String \\
            -a clairs-to \\
            -o tmp.vcf.gz
    fi

    clairs_to_fixes.bash tmp.vcf.gz ${output}

    get_nextflow_log.bash clairs-to.log
    """
}
// TODO: can use panel of normals
