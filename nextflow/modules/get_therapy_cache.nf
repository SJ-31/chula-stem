process GET_THERAPY_CACHE {
    ext version: "1"

    input:
    tuple val(meta), path(files)

    output:
    //
    path ("therapy_cache/${civic}"), emit: civic
    path ("therapy_cache/${pandrugs2}"), emit: pandrugs2
    path not_in
    tuple val(meta), path("*.log"), emit: log

    script:
    //
    civic = "civic_cache.json"
    pandrugs2 = "pandrugs2_cache.json"
    not_in = "not_in_databases.txt"
    check = file("${meta.out}/therapy_cache")
    if (check) {
        """
        cp -r ${check} .
        ln -sr ${meta.log}/therapy_cache.log .
        """
    }
    else {
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
