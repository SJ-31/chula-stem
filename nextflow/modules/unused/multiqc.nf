process MULTIQC {
    ext version: "1.25.1"

    input:
    tuple val(meta), path(files)

    output:

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
    path (out), emit: html
    path (data), emit: data
    path ("multiqc.log"), emit: logfile

    script:
    //
    out = Utl.getName(meta, "MultiQC", "html")
    data = Utl.getName(meta, "MultiQC_data")
    """
    multiqc . --config "${params.configdir}/multiqc_config.yaml" \\
        --filename ${out}

    get_nextflow_log.bash multiqc.log
    """
}
