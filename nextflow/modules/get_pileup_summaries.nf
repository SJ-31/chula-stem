process GET_PILEUP_SUMMARIES {
    ext version: params.gatk_version

    label "big_mem"
    publishDir "$meta.out", mode: "copy", saveAs: { x -> x ==~ /.*\.log/ ? null : x }
    publishDir "$meta.log", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(bam), path(indices, arity: "1..*")
    val(variants) // Commmon germline variant sites, with only biallelic SNPs
    val(module_number)

    //

    output:
    tuple val(meta), path(output), emit: pileup
    path("*.log")
    //

    script:
    output = "${module_number}-${meta.filename}-Pileup.table"
    check = file("${meta.out}/${output}")
    args = task.ext.args.join(" ")
    if (check.exists()) {
        """
        ln -sr ${check} .
        ln -sr ${meta.log}/get_pileup_summaries.log .
        """
    } else {
        """
        gatk GetPileupSummaries \\
            -I ${bam} \\
            -L ${variants} \\
            -V ${variants} \\
            -O ${output}

        get_nextflow_log.bash get_pileup_summaries.log
        """
    }
    //
}
