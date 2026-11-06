process SAREK_MANIFEST {

    input:
    val(sarek_config)
    // TODO: [2026-11-06 Fri] two types of outputs: the paths to trimmed
    // fastq files
    // and recalibrated bam, labelled with tumor, normal
    // and the sarek outputs, labelled by group: germline, tumor_only,
    // and paired

    output:
    path("tool_outs.csv"), emit: tool_outs
    path("preprocessed.csv"), emit: preprocessed

    script:
    """
    make_sarek_manifest -c ${sarek_config} \\
       -t tool_outs.csv \\
       -p preprocessed.csv
    """
}
