process PURECN_NORMALDB {
    ext version: "2.16.1"

    input:
    tuple val(meta), path(normals)
    val pon

    output:

    //
    path (output), emit: db
    path (mapping_bias), emit: mapping_bias, optional: true
    path (mapping_bias_hq), optional: true
    path png
    path low_cov
    tuple val(meta), path("*.log"), emit: log

    script:
    //
    output = "normalDB_${params.genome_build}.rds"
    png = "interval_weights_${params.genome_build}.png"
    low_cov = "low_coverage_targets_${params.genome_build}.bed"
    mapping_bias = "mapping_bias_hg38.rds"
    mapping_bias_hq = "mapping_bias_hq_sites_hg38.bed"

    args = task.ext.args.join(" ")
    pon_flag = pon ? " --normal-panel ${pon}" : ""

    """
    ls *coverage*.txt.gz > files.list

    Rscript \$PURECN/NormalDB.R \\
        ${args} \\
        --out-dir . \\
        --coverage-files files.list \\
        ${pon_flag} \\
        --genome ${params.genome_build}

    get_nextflow_log.bash purecn_normaldb.log
    """
    //
}
