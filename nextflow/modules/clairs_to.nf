process CLAIRS_TO {
    ext version: "0.3.1"

    publishDir "${meta.out}", mode:"copy", saveAs: params.saveFn
    publishDir "${meta.log}", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(tumor), path(index)
    val(reference)
    val(target_intervals)
    val(module_number)
    //

    output:
    tuple val(meta), path(check), emit: variants
    path("*.log")
    //

    script:
    output = Utils.getName(module_number, meta, "ClairS-TO", "vcf.gz")
    check = file("${meta.out}/${output}")
    target_flag = target_intervals != "" ? " --bed_fn ${target_intervals} " : ""
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
            --threads ${task.cpus} \\
            --platform ${task.ext.platform} \\
            ${target_flag} \\
            --output_dir .

        bcftools concat snv.vcf.gz indel.vcf.gz -a | \\
            vcf_info_add_tag.bash -n ${params.source_name} \\
                -d "$params.source_description" \\
                -b '.' \\
                -t String \\
                -a clairs-to \\
                -o ${output}

        get_nextflow_log.bash clairs-to.log
        """
    }
    //
}
// TODO: can use panel of normals
