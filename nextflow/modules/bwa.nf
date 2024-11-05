process BWA {
    publishDir "$meta.out", mode: "copy"
    publishDir "$meta.log", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), val(reads)
    val(reference)
    val(module_number)
    //

    output:
    path("*.log")
    //

    script:
    """
    bwa-mem2 index $reference
    bwa-mem2 mem \
        -o "${module_number}-${meta.id}.sam" \
        $reference \
        ${reads[0]} ${reads[1]}

    cp .command.out bwa.log
    """
    //
}
