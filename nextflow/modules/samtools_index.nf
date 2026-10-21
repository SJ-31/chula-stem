process SAMTOOLS_INDEX {
    ext version: "1.21"

    
    

    input:
    tuple val(meta), path(bam)
    //

    output:
    tuple val(meta), path(output), emit: index
    tuple val(meta), path("*.log"), emit: log
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
