process SAREK_MANIFEST {

    input:
    val(sarek_config)

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
