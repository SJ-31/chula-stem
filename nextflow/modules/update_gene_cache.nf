nextflow.enable.types = true

process UPDATE_GENE_CACHE {
    input:
    vcfs: Bag<Path>
    cache_file: String
    open_targets: String?
    gene_column: String
    filter_confident: Boolean
    

    output:
    record(cache: cache_file, failed_lookups: file("failed_lookups.txt", optional: true))
    
    script:
    input_str = vcfs.collect { v -> "-i ${v}" }.join(" ")
    filter_str = filter_confident ? "--filter-confident" : ""
    open_targets_str = open_targets ? "export OPEN_TARGETS_DUMP=${open_targets}" : ""
    """
    ${open_targets_str}
    update_gene_cache ${input_str} ${filter_str} \\
        -e ${cache_file} \\
        -c ${gene_column} \\
        -d open_targets -d civic \\
        --not-in-db failed_lookups.txt
    """
        
}
