process GET_PILEUP_SUMMARIES {
    ext version: params.gatk_version

    label "big_mem"
    
    

    input:
    tuple val(meta), path(bam), path(indices, arity: "1..*")
    val(variants) // Commmon germline variant sites, with only biallelic SNPs
    

    //

    output:
    tuple val(meta), path(output), emit: pileup
    tuple val(meta), path("*.log"), emit: log
    //

    script:
    output = Utl.getName(meta, "Pileup", "table")
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
            --java-options "-Xmx8G" \\
            -I ${bam} \\
            -L ${variants} \\
            -V ${variants} \\
            -O ${output}

        get_nextflow_log.bash get_pileup_summaries.log
        """
    }
    //
}
