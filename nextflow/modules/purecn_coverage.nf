process PURECN_COVERAGE {
    ext version: ""

    input:
    tuple val(meta), path(bam_or_cov), path(index)
    val bait_intervals
    val normalize_only // PureCN supports using third-party tools for coverage calculation. If `normalize_only` is true, then the precomputed coverage files are GC-normalized

    output:

    //
    tuple val(meta), path(o1), emit: loess
    tuple val(meta), path(o2), emit: cov
    tuple val(meta), path("*.log"), emit: log

    script:
    //
    o1 = Utl.getName(meta, "coverage_loess", "txt.gz")
    o2 = Utl.getName(meta, "coverage", "txt.gz")
    sample_flag = normalize_only ? " --coverage " : " --bam "
    args = task.ext.args.join(" ")

    """
    Rscript \$PURECN/Coverage.R \\
        ${args} \\
        --out-dir . \\
        ${sample_flag} ${bam_or_cov} \\
        --intervals ${bait_intervals} \\
        --cores ${task.cpus}

    mv *_loess.txt.gz ${o1}
    mv *_coverage.txt.gz ${o2}

    get_nextflow_log.bash purecn_coverage.log
    """
}
