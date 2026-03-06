process PURECN_CALL {
    ext version: "2.16.1"

    publishDir "${meta.out}", mode:"copy", saveAs: params.saveFn
    publishDir "${meta.log}", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(coverage), path(mutect2_vcf), path(mutect2_stats)
    path(normaldb)
    path(bait_intervals)
    path(mapping_bias)
    val(snp_blacklist)
    val(module_number)
    //

    output:
    // TODO: still unsure about the output
    tuple val(meta), path(o1), emit: loess
    path("*.log")
    //

    script:
    // o1 = Utl.getName(module_number, meta, "PureCN_coverage_loess", ".txt.gz")
    // c1 = file("${meta.out}/${o1}")
    blacklist_flag = snp_blacklist ? " --snp-blacklist ${snp_blacklist}" : ""
    args = task.ext.args.join(" ")
    if (c1.exists()) {
        """
        ln -sr ${c1} .
        ln -sr ${c2} .
        ln -sr ${meta.log}/purecn_coverage.log .
        """
    } else {
        """
        Rscript "${params.purecn_extdata}/PureCN.R" \\
            ${args} \\
            ${blacklist_flag} \\
            --out . \\
            --tumor "${coverage}" \\
            --vcf "${mutect2_vcf}" \\
            --stats-file "${mutect2_stats}" \\
            --normaldb "${normaldb}" \\
            --intervals "${bait_intervals}" \\
            --genome ${params.genome_build}


        get_nextflow_log.bash purecn_call.log
        """
    }
    //
}
