process MOSDEPTH {
    ext version: "0.3.9"
 

    publishDir "$meta.out", mode: 'copy'
    publishDir "$meta.log", mode: 'copy', pattern: "*.log"

    input:
    tuple val(meta), path(bam)
    val(target_regions)
    val(module_number)

    output:
    tuple val(meta), path("${prefix}*")
    //

    script:
    prefix = "${module_number}-${meta.id}_Mosdepth"
    check = file("${meta.out}/${prefix}.mosdepth.global.dist.txt")
    if (check.exists()) {
        """
        cp "${meta.out}/${prefix}*" .
        cp "${meta.log}/mosdepth.log" .
        """
    } else if (target_regions) {
        """
        mosdepth $prefix $bam \\
            --by $target_regions
        cp .command.out mosdepth.log
        """
    } else {
        """
        mosdepth $prefix $bam
        cp .command.out mosdepth.log
        """
    }

    //
}
