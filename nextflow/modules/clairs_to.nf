process CLAIRS_TO {
    ext version: "0.3.1"

    publishDir "${meta.out}", mode:"copy", saveAs: params.saveFn
    publishDir "${meta.log}", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(tumor), path(index), path(previous_variants)
    val(reference)
    val(target_intervals)
    val(module_number)
    //

    output:
    tuple val(meta), path(output), emit: variants
    path("*.log")
    //

    script:
    output = Utl.getName(module_number, meta, "ClairS-TO", "vcf.gz")
    check = file("${meta.out}/${output}")
    target_flag = target_intervals != "" ? " --bed_fn ${target_intervals} " : ""
    prev_flag = file(previous_variants).size() != 0 ? "--hybrid_mode_vcf_fn ${previous_variants}" : ""
    args = task.ext.args.join(" ")
    if (check.exists()) {
        """
        ln -sr ${check} .
        ln -sr ${meta.log}/clairs-to.log .
        """
    } else {
        """
        run_clairs_to ${args} \\
            --tumor_bam_fn ${tumor} \\
            --ref_fn ${reference} \\
            --sample_name ${meta.RGSM_tumor} \\
            --threads ${task.cpus} \\
            --platform ${task.ext.platform} \\
            ${target_flag} \\
            ${prev_flag} \\
            --output_dir .

        if [[ -e snv.vcf.gz && -e indel.vcf.gz ]]; then
            bcftools concat snv.vcf.gz indel.vcf.gz -a | \\
                vcf_info_add_tag.bash -n ${params.source_name} \\
                    -d "$params.source_description" \\
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
                -d "$params.source_description" \\
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
    //
}
// TODO: can use panel of normals
