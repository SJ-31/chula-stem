process PURECN_BAIT_INTERVALS {
    ext version: "2.16.1"

    publishDir "${meta.out}", mode:"copy", saveAs: params.saveFn
    publishDir "${meta.log}", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), val(genome)
    val(baits)
    val(mappability)
    val(module_number)
    //

    output:
    path(o1), emit: baits
    path(o2), emit: baits_optimized
    path("*.log")
    //

    script:
    o1 = "baits_intervals.txt"
    o2 = "baits_optimized.bed"
    c1 = file("${meta.out}/${o1}")
    c2 = file("${meta.out}/${o2}")
    args = task.ext.args.join(" ")
    if (c1.exists()) {
        """
        ln -sr ${c1} .
        ln -sr ${c2} .
        ln -sr ${meta.log}/purecn_bait_intervals.log .
        """
    } else {
        """
        Rscript \$PURECN/IntervalFile.R \\
        ${args} \\
        --in-file "${baits}" \\
        --fasta "${genome}" \\
        --genome ${params.genome_build} \\
        --out-file "${o1}" \\
        --export "${o2}" \\
        --mappability "${mappability}"

        get_nextflow_log.bash purecn_bait_intervals.log
        """
    }
    //
}
