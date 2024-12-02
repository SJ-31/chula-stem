process MULTIQC {
    ext version: "1.25.1"

    publishDir "$meta.out", mode: "copy", saveAs: params.saveFn
    publishDir "$meta.log", mode: 'copy', pattern: "*.log"

    input:
    tuple val(meta), path(files)
    val(module_number)
    // Will aggregate metrics from
    // - fastp
    //      from json file
    // - Picard (AlignmentSummaryMetrics,
    //                  CollectHsMetrics|CollectWgsMetrics|CollectRnaSeqMetrics)
    //                  from
    // - mosdepth (coverage of reference, sequencing depth)
    //    Ideally "region.dist.txt", but "global.dist.txt" if region is unavailable
    // - bcftools stats
    //    from plot.py
    // - VEP
    //  from html summary file

    output:
    path("multiqc")
    //

    script:
    check = file("${meta.out}/multiqc")
    if (check.exists()) {
        """
        ln -sr -r $check .
        ln -sr "${meta.log}/multiqc.log" .
        """
    } else {
        """
        multiqc . --config "$params.configdir/multiqc_config.yaml"

        get_nextflow_log.bash multiqc.log
        """
    }

    //
}
