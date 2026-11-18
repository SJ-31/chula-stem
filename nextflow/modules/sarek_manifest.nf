nextflow.enable.types = true

process SAREK_MANIFEST {
    input:
    sarek_config: Path

    output:
    record(tool_outs: file("tool_outs.csv"), preprocessed: file("preprocessed.csv"))

    script:
    """
    make_sarek_manifest -c ${sarek_config} \\
       -t tool_outs.csv \\
       -p preprocessed.csv
    """
}
