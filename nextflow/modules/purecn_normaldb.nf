nextflow.enable.types = true

process PURECN_NORMALDB {
    ext version: "2.16.1"

    input:
    normals: Bag<Path>
    pon: String
    genome_build: String

    output:
    record(
        db: file(output),
        mapping_bias: file(mapping_bias),
        mapping_bias_hq: file(mapping_bias_hq),
        png: file(png),
        low_cov: low_cov,
    )

    topic:
    record(process: "purecn_normaldb", file: file("purecn_normaldb.log")) >> "logs"

    script:
    output = "normalDB_${genome_build}.rds"
    png = "interval_weights_${genome_build}.png"
    low_cov = "low_coverage_targets_${genome_build}.bed"
    mapping_bias = "mapping_bias_hg38.rds"
    mapping_bias_hq = "mapping_bias_hq_sites_hg38.bed"

    args = task.ext.args ?: ""
    pon_flag = pon ? " --normal-panel ${pon}" : ""

    """
    ls *coverage*.txt.gz > files.list

    Rscript \$PURECN/NormalDB.R \\
        ${args} \\
        --out-dir . \\
        --coverage-files files.list \\
        ${pon_flag} \\
        --genome ${genome_build}

    get_nextflow_log.bash purecn_normaldb.log
    """
    //
}
