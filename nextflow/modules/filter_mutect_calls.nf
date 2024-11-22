process FILTER_MUTECT_CALLS {
    ext version: params.gatk_version

    publishDir "$meta.out", mode: "copy", saveAs: { x -> x ==~ /.*\.log/ ? null : x }
    publishDir "$meta.log", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(vcf), path(ro_model), path(contamination), path(segmentation), path(stats)
    val(reference)
    val(module_number)
    //

    output:
    tuple val(meta), path(output), emit: filtered
    path("*.log")
    //

    script:
    output = "${module_number}-${meta.filename}-Mutect2_filtered.vcf.gz"
    check = file("${meta.out}/${output}")
    args = task.ext.args.join(" ")
    if (check.exists()) {
        """
        ln -sr ${check} .
        ln -sr ${meta.log}/filter_mutect_calls.log .
        """
    } else {
        """
        gatk IndexFeatureFile -I ${vcf}

        gatk FilterMutectCalls -V ${vcf} \\
            --ob-priors ${ro_model} \\
            --tumor-segmentation ${segmentation} \\
            --contamination-table ${contamination} \\
            --reference ${reference} \\
            -O ${output}

        get_nextflow_log.bash filter_mutect_calls.log
        """
    }
    //
}
