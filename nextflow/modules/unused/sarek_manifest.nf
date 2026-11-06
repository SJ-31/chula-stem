process SAREK_MANIFEST {
    ext version: "1.0.0"

    input:
    val(sarek_config)

    output:
    path("germline.csv"), emit: germline
    path("tumor_only.csv"), emit: tumor_only
    path("paired.csv"), emit: paired

    script:
    """
    make_sarek_manifest -c ${sarek_config} \\
       -g germline.csv \\
       -t tumor_only.csv \\
       -p paired.csv
    """
}
