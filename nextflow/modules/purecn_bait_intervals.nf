process PURECN_BAIT_INTERVALS {
    ext version: "2.16.1"

    input:
    tuple val(meta), val(genome)
    val baits
    val mappability

    output:

    //
    path (o1), emit: baits
    path (o2), emit: baits_optimized
    tuple val(meta), path("*.log"), emit: logfile

    script:
    //
    o1 = "baits_intervals.txt"
    o2 = "baits_optimized.bed"
    args = task.ext.args.join(" ")
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
