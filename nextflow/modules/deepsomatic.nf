process DEEPSOMATIC {
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
    output = Utils.getName(module_number, meta, "Deepsomatic", "vcf.gz")
    check = file("${meta.out}/${output}")
    target_flag = target_intervals != "" ? " --regions ${target_intervals} " : ""
    args = task.ext.args.join(" ")
    if (check.exists()) {
        """
        ln -sr ${check} .
        ln -sr ${meta.log}/deepsomatic.log .
        """
    } else {
        """
        run_deepsomatic ${args} \\
            --model_type ${task.ext.model} \\
            --reads_normal ${normal} \\
            --reads_tumor ${tumor} \\
            --ref ${reference} \\
            ${target_flag} \\
            --sample_name_normal ${meta.RGSM_normal} \\
            --sample_name_tumor ${meta.RGSM_tumor} \\
            --process_somatic true \\
            --output_vcf temp.vcf.gz

        echo -e "${meta.RGSM_tumor}\t${meta.RGSM_normal}" > rename.txt

        bcftools filter -i 'FILTER~\"GERMLINE\"' temp.vcf.gz -O z | \\
            bcftools reheader -s rename.txt > germline.vcf.gz

        bcftools index germline.vcf.gz

        bcftools filter -i 'FILTER!~\"GERMLINE\"' temp.vcf.gz -O z > somatic.vcf.gz

        bcftools index somatic.vcf.gz

        bcftools merge germline.vcf.gz somatic.vcf.gz | \\
            bcftools view -s "${meta.RGSM_normal},${meta.RGSM_tumor}" | \\
            vcf_info_add_tag.bash -n ${params.source_name} \\
                -d "$params.source_description" \\
                -b '.' \\
                -t String \\
                -a deepsomatic \\
                -o ${output}

        get_nextflow_log.bash deepsomatic.log
        """
    }
    //
}
// TODO: can be used with panel of normals
