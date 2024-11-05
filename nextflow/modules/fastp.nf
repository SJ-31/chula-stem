process FASTP {
    publishDir "$outdir/$patient", mode: "copy"
    publishDir "$logdir/$patient", mode: "copy", pattern: "*.log"

    input:
    tuple val(patient), path(reads)
    val(outdir)
    val(logdir)
    //

    output:
    path("*.log")
    //

    script:
    """


    cp .command.out .log
    """
    //
}
