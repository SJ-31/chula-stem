process PURECN_COVERAGE {
    ext version: ""

    publishDir "${meta.out}", mode:"copy", saveAs: params.saveFn
    publishDir "${meta.log}", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(bam_or_cov), path(index)
    val(bait_intervals)
    val(normalize_only) // PureCN supports using third-party tools for coverage calculation. If `normalize_only` is true, then the precomputed coverage files are GC-normalized
    val(module_number)
    //

    output:
    tuple val(meta), path(o1), emit: loess
    tuple val(meta), path(o2), emit: cov
    path("*.log")
    //

    script:
    o1 = Utl.getName(module_number, meta, "coverage_loess", "txt.gz")
    o2 = Utl.getName(module_number, meta, "coverage", "txt.gz")
    c1 = file("${meta.out}/${o1}")
    c2 = file("${meta.out}/${o2}")
    sample_flag = normalize_only ? " --coverage " : " --bam "
    args = task.ext.args.join(" ")
    if (c1.exists()) {
        """
        ln -sr ${c1} .
        ln -sr ${c2} .
        ln -sr ${meta.log}/purecn_coverage.log .
        """
    } else {
        """
        Rscript \$PURECN/Coverage.R \\
            ${args} \\
            --out-dir . \\
            ${sample_flag} ${bam_or_cov} \\
            --intervals ${bait_intervals} \\
            --cores ${task.cpus} > tmp

        mv *_loess.txt.gz ${o1}
        mv *_coverage.txt.gz ${o2}
        
        get_nextflow_log.bash purecn_coverage.log
        cat tmp >> purecn_coverage.log
        """
    }
    //
}
