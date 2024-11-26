process SAMTOOLS_INDEX {
    ext version: "1.21"

    publishDir "$meta.out", mode: "copy", saveAs: params.saveFn
    publishDir "$meta.log", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(bam)
    //

    output:
    tuple val(meta), path(output), emit: index
    path("*.log")
    //

    script:
    output = "${bam}.bai"
    check = file("${meta.out}/$output")
    args = task.ext.args.join(" ")
    if (check.exists()) {
        """
        ln -sr $check .
        ln -sr ${meta.log}/samtools_index.log .
        """
    } else {
        """
        samtools index $bam
        get_nextflow_log.bash samtools_index.log
        """
    }
    //
}
