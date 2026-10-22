process MOSDEPTH {
    ext version: "0.3.9"

    input:
    tuple val(meta), path(bam), path(index)
    val target_regions

    output:
    tuple val(meta), path("${prefix}*"), emit: all
    tuple val(meta), path(out), emit: dist
    tuple val(meta), path("*log"), emit: logfile

    script:
    //
    prefix = Utl.getName(meta, "Mosdepth")
    if (target_regions) {
        out = "${prefix}.mosdepth.region.dist.txt"
    }
    else {
        out = "${prefix}.mosdepth.global.dist.txt"
    }
    if (target_regions) {
        """
        mosdepth ${prefix} ${bam} \\
            --by ${target_regions}
        get_nextflow_log.bash mosdepth.log
        """
    }
    else {
        """
        mosdepth ${prefix} ${bam}
        get_nextflow_log.bash mosdepth.log
        """
    }

    //
}
