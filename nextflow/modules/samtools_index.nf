process SAMTOOLS_INDEX {
    ext version: "1.21"

    publishDir "$meta.out", mode: "copy"
    publishDir "$meta.log", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(bam)
    val(meta)
    //

    output:
    tuple val(meta), path(output), emit: index
    path("*.log")
    //

    script:
    output = "${bam.baseName}.bai"
    check = file("${meta.out}/$output")
    args = task.ext.args.join(" ")
    if (check.exists()) {
        """
        cp $check .
        cp ${meta.log}/samtools_index .
        """
    } else {
        """
        samtools index $bam
        cp .command.out samtools_index.log
        """
    }
    //
}
