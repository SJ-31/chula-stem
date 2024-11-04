process FASTP {
    publishDir "$outdir", mode: "copy"
    publishDir "$params.logdir", mode: "copy", pattern: "*.log"

    input:
    tuple val(patient), path(reads)
    val(outdir)
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
