process CLAIRS {
    ext version: "0.4.1"

    publishDir "${meta.out}", mode:"copy", saveAs: params.saveFn
    publishDir "${meta.log}", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(normal), path(tumor), path(indices, arity: "2")
    val(reference)
    val(target_intervals)
    val(module_number)
    //

    output:
    tuple val(meta), path(output)
    path("*.log")
    //

    script:
    output = Utils.getName(module_number, meta, "ClairS", "vcf.gz")
    check = file("${meta.out}/${output}")
    target_flag = target_intervals != "" ? " --bed_fn ${target_intervals} " : ""
    args = task.ext.args.join(" ")
    if (check.exists()) {
        """
        ln -sr ${check} .
        ln -sr ${meta.log}/clairs.log .
        """
    } else {
        """
        run_clairs ${args} \\
            --tumor_bam_fn ${tumor} \\
            --normal_bam_fn ${normal} \\
            --ref_fn ${reference} \\
            --threads ${task.cpus} \\
            --platform ${task.ext.platform} \\
            --sample_name ${meta.RGSM_tumor} \\
            ${target_flag} \\
            --output_dir .

        echo -e "${meta.RGSM_tumor}\t${meta.RGSM_normal}" > rename.txt
        bcftools filter -i 'FILTER~\"Germline\"' temp.vcf.gz -O z | \\
            bcftools reheader -s rename.txt > germline.vcf.gz

        bcftools index germline.vcf.gz
        bcftools filter -i 'FILTER!~\"Germline\"' temp.vcf.gz -O z > somatic.vcf.gz
        bcftools index somatic.vcf.gz

        bcftools merge germline.vcf.gz somatic.vcf.gz | \\
            bcftools view -s "${meta.RGSM_normal},${meta.RGSM_tumor}" | \\
            vcf_info_add_tag.bash -n ${params.source_name} \\
                -d "$params.source_description" \\
                -b '.' \\
                -t String \\
                -a clairs \\
                -o ${output}

        get_nextflow_log.bash clairs.log
        """
    }
    //
}
