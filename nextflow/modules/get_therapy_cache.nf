process GET_THERAPY_CACHE {
    ext version: "1"

    publishDir "${meta.out}", mode:"copy", saveAs: params.saveFn
    publishDir "${meta.log}", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(files)
    //

    output:
    path("therapy_cache/${civic}"), emit: civic
    path("therapy_cache/${pandrugs2}"), emit: pandrugs2
    path(not_in)
    path("*.log")
    //

    script:
    civic = "civic_cache.json"
    pandrugs2 = "pandrugs2_cache.json"
    not_in = "not_in_databases.txt"
    check = file("${meta.out}/therapy_cache")
    if (check) {
        """
        cp -r ${check} .
        ln -sr ${meta.log}/therapy_cache.log .
        """
    } else {
        """
        mkdir therapy_cache

        get_therapy_cache -i . \\
            --filter-confident \\
            --not-in-db therapy_cache/${not_in} \\
            --civic-cache therapy_cache/${civic} \\
            --pandrugs2-cache therapy_cache/${pandrugs2} \\
            --gene-col SYMBOL

        get_nextflow_log.bash therapy_cache.log
        """
    }
    //
}
