nextflow.enable.types = true

process PURECN_BAIT_INTERVALS {
    ext version: "2.16.1"

    input:
    genome: Path
    baits: Path
    mappability: Path

    output:
    record(baits: file(o1),
           baits_optimized: file(o2)
    )

    topic:
    record(process: "purecn_bait_intervals", file: file("purecn_bait_intervals.log")) >> "logs"

    script:
    //
    o1 = "baits_intervals.txt"
    o2 = "baits_optimized.bed"
    args = task.ext.args ?: ""
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
