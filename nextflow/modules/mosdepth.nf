process MOSDEPTH {
    ext version: "0.3.9"
 

    publishDir "$meta.out", mode: "copy", saveAs: params.saveFn
    publishDir "$meta.log", mode: 'copy', pattern: "*.log"

    input:
    tuple val(meta), path(bam)
    val(target_regions)
    val(module_number)

    output:
    tuple val(meta), path("${prefix}*")
    //

    script:
    prefix = "${module_number}-${meta.filename}-Mosdepth"
    check = file("${meta.out}/${prefix}.mosdepth.global.dist.txt")
    if (check.exists()) {
        """
        ln -sr "${meta.out}/${prefix}*" .
        ln -sr "${meta.log}/mosdepth.log" .
        """
    } else if (target_regions) {
        """
        mosdepth $prefix $bam \\
            --by $target_regions
        get_nextflow_log.bash mosdepth.log
        """
    } else {
        """
        mosdepth $prefix $bam
        get_nextflow_log.bash mosdepth.log
        """
    }

    //
}
