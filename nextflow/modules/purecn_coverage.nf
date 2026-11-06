nextflow.enable.types = true

process PURECN_COVERAGE {
    ext version: ""

    input:
    record(id: String,
           bam_or_cov: Path,
           index: Path
    )
    bait_intervals: Path
    normalize_only: Boolean
    // PureCN supports using third-party tools for coverage calculation. If `normalize_only` is true, then the precomputed coverage files are GC-normalized

    output:
    record(id: id,
           loess: file(o1),
           coverage: file(o2)
    )

    topic:
    record(id: id, process: "purecn_coverage", file: "purecn_coverage.log") >> "logs"
    
    script:
    o1 = "${id}-coverage_loess.txt.gz"
    o2 = "${id}-coverage.txt.gz"
    sample_flag = normalize_only ? " --coverage " : " --bam "
    args = task.ext.args ?: ""

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
