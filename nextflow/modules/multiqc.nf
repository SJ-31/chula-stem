process MULTIQC {
    ext version: "1.25.1"

    publishDir "$meta.out", mode: "copy", saveAs: { x -> x ==~ /.*\.log/ ? null : x }
    publishDir "$meta.log", mode: 'copy', pattern: "*.log"

    input:
    tuple val(meta), path(metrics)
    val(module_number)
    // Will aggregate metrics from
    // - fastp
    // - Picard & GATK (AlignmentSummaryMetrics, MarkDuplicates,
    //                  CollectHsMetrics|CollectWgsMetrics|CollectRnaSeqMetrics)
    // - mosdepth (coverage of reference, sequencing depth)

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
