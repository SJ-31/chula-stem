process GET_THERAPY_CACHE {
    ext version: "1"

    publishDir "${meta.out}", mode:"copy", saveAs: params.saveFn
    publishDir "${meta.log}", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(files)
    //

    output:
    path(civic), emit: civic
    path(pandrugs2), emit: pandrugs2
    path("*.log")
    //

    script:
    civic = "civic_cache.json"
    pandrugs2 = "pandrugs2_cache.json"
    not_in = "not_in_databases.txt"
    check1 = file("${meta.out}/${civic}")
    check2 = file("${meta.out}/${pandrugs2}")
    if (check1.exists() && check2.exists()) {
        """
        ln -sr ${check1} .
        ln -sr ${check2} .
        ln -sr ${meta.out}/${not_in} .
        ln -sr ${meta.log}/therapy_cache.log .
        """
    } else {
        """
        get_therapy_cache -i . \\
            --filter_confident \\
            --not_in_db ${not_in} \\
            --civic_cache ${civic} \\
            --pandrugs2_cache ${pandrugs2}

        get_nextflow_log.bash therapy_cache.log
        """
    }
    //
}
