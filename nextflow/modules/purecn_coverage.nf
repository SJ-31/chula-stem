process PURECN_COVERAGE {
    ext version: ""

    publishDir "${meta.out}", mode:"copy", saveAs: params.saveFn
    publishDir "${meta.log}", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(bam_or_cov)
    val(bait_intervals)
    val(normalize_only) // PureCN supports using third-party tools for coverage calculation. If `normalize_only` is true, then the precomputed coverage files are GC-normalized
    val(module_number)
    //

    output:
    path("*.log")
    //

    script:
    output = Utl.getName(module_number, meta, "PureCN_loess", ".txt.gz")
    check = file("${meta.out}/${output}")
    sample_flag = normalize_only ? " --coverage " : " --bam "
    args = task.ext.args.join(" ")
    if (check.exists()) {
        """
        ln -sr ${check} .
        ln -sr ${meta.log}/purecn_coverage.log .
        """
    } else {
        """
        Rscript ${params.purecn_extdata}/Coverage.R \\
            ${args} \\
            --out-dir . \\
            ${sample_flag} ${bam_or_cov} \\
            --intervals ${bait_intervals} \\
            --cores ${task.cpus}

        get_nextflow_log.bash purecn_coverage.log
        """
    }
    //
}
