process OCTOPUS {
    ext version: "0.7.4"

    label "big_mem"
    publishDir "$meta.out", mode:"copy", saveAs: params.saveFn
    publishDir "$meta.log", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(normal), path(tumor), path(indices, arity: "2"), path(previous_variants)
    // Previously called variants will be added to
    //      octopus' list of candidate variants (does not guarantee that they will be called)
    val(reference) // Must have ".fasta" extension
    val(target_intervals) // (bgzipped) BED file of target intervals
    val(tumor_only)
    val(module_number)
    //

    output:
    tuple val(meta), path(somatic), emit: variants
    path(output)
    path("*.log")
    //

    script:
    output = Utl.getName(module_number, meta, "Octopus_all", "vcf.gz")
    somatic = Utl.getName(module_number, meta, "Octopus", "vcf.gz")
    check1 = file("${meta.out}/${output}")
    prev_flag = previous_variants.size() != 0 ? "--source-candidates ${previous_variants}" : ""
    if (!tumor_only) {
        bam_flag = "-I ${normal} ${tumor} --normal-sample ${meta.RGSM_normal}"
    } else {
        bam_flag = "-I ${tumor} -C cancer"
    }
    targets_flag = target_intervals ? "--regions-file ${target_intervals} " : ""
    args = task.ext.args.join(" ")
    if (check1.exists()) {
        """
        ln -sr ${check1} .
        ln -sr ${meta.log}/octopus.log .

        bcftools filter -i "INFO/SOMATIC == 1" ${output} -O z > ${somatic}
        """
    } else {
        """
        if [[ -e "${previous_variants}" && -n "${prev_flag}" ]]; then
            bcftools index ${previous_variants}
        fi

        octopus --reference ${reference} \\
            ${bam_flag} \\
            ${prev_flag} \\
            ${targets_flag} \\
            ${args} \\
            --threads ${task.cpus} \\
            --output tmp.vcf

        if [[ "${tumor_only}" == "false" ]]; then
            bcftools view -s "${meta.RGSM_normal},${meta.RGSM_tumor}" tmp.vcf | \\
                vcf_info_add_tag.bash -n ${params.source_name} \\
                    -d "$params.source_description" \\
                    -b '.' \\
                    -t String \\
                    -a octopus \\
                    -o ${output}
        else
            vcf_info_add_tag.bash -n ${params.source_name} \\
                -d "$params.source_description" \\
                -b '.' \\
                -t String \\
                -a octopus \\
                -i tmp.vcf \\
                -o ${output}
        fi
        bcftools filter -i "INFO/SOMATIC == 1" ${output} -O z > ${somatic}

        get_nextflow_log.bash octopus.log
        """
    }
    //
}

