process GET_THERAPY_CACHE {
    ext version: "1"

    input:
    path small_variants
    path structural_variants
    path cache_path

    output:
    //
    path cache_path, emit: cache
    path not_in, emit: failed
    path("*.log"), emit: logfile

    script:
    //
    not_in = "not_in_databases.txt"
    """
    
    
    get_therapy_cache -i . \\
         --filter-confident \\
         --not-in-db therapy_cache/${not_in} \\
         --cache ${cache_path} \\
         --gene-col SYMBOL

    get_nextflow_log.bash therapy_cache.log
    """
}
