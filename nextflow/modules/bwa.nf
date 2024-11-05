process BWA {
    // bwa-mem2 version 2.2.1
    publishDir "$meta.out", mode: "copy"
    publishDir "$meta.log", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), val(reads)
    val(reference)
    val(module_number)
    //

    output:
    tuple val(meta), path(out), emit: mapped
    path("*.log")
    //

    script:
    out = "${module_number}-${meta.id}.sam"
    """
    bwa-mem2 index $reference
    bwa-mem2 mem \
        -o $out \
        -v 3 \
        $reference \
        ${reads[0]} ${reads[1]}

    cp .command.out bwa.log
    """
    //
}
