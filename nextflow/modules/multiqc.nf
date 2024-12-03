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
    path(out)
    path(data)
    path("multiqc.log")
    //

    script:
    out = Utils.getName(module_number, meta, "MultiQC", "html")
    data = Utils.getName(module_number, meta, "MultiQC_data")
    check = file("${meta.out}/${out}")
    check2 = file("${meta.out}/${data}")
    if (check.exists()) {
        """
        cp -r ${check2} .
        ln -sr ${check} .
        ln -sr "${meta.log}/multiqc.log" .
        """
    } else {
        """
        multiqc . --config "$params.configdir/multiqc_config.yaml" \\
            --filename ${out}

        get_nextflow_log.bash multiqc.log
        """
    }

    //
}
