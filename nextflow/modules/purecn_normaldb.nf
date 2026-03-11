process PURECN_NORMALDB {
    ext version: "2.16.1"

    publishDir "${meta.out}", mode:"copy", saveAs: params.saveFn
    publishDir "${meta.log}", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(normals)
    val(pon)
    val(module_number)
    //

    output:
    path(output), emit: db
    path(mapping_bias), emit: mapping_bias, optional: true
    path(mapping_bias_hq), optional: true
    path(png)
    path(low_cov)
    path("*.log")
    //

    script:
    output = "normalDB_${params.genome_build}.rds"
    png = "interval_weights_${params.genome_build}.png"
    low_cov = "low_coverage_targets_${params.genome_build}.bed"
    mapping_bias = "mapping_bias_hg38.rds"
    mapping_bias_hq = "mapping_bias_hq_sites_hg38.bed"

    c1 = file("${meta.out}/${output}")
    c2 = file("${meta.out}/${png}")
    c3 = file("${meta.out}/${low_cov}")
    c4 = file("${meta.out}/${mapping_bias}")
    c5 = file("${meta.out}/${mapping_bias_hq}")

    args = task.ext.args.join(" ")
    pon_flag = pon ? " --normal-panel ${pon}" : ""
    if (c1.exists()) {
        """
        ln -sr ${c1} .
        ln -sr ${c2} .
        ln -sr ${c3} .
        ln -sr ${meta.log}/purecn_normaldb.log .

        if [[ -e ${c4} ]]; then
            ln -sr ${c4} .
        fi
        if [[ -e ${c5} ]]; then
            ln -sr ${c5} .
        fi
        """
    } else {
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
    }
    //
}
