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
    path(output)
    path("*.log")
    //

    script:
    output = "${module_number}-normalDB_${params.genome_build}.rds"
    check = file("${meta.out}/${output}")
    args = task.ext.args.join(" ")
    pon_flag = !pon ? " --normal-panel ${pon}" : ""
    if (check.exists()) {
        """
        ln -sr ${check} .
        ln -sr ${meta.log}/ .
        """
    } else {
        """
        ls *coverage.txt.gz > files.list

        Rscript "${params.purecn_extdata}/NormalDB.R"
            ${args} \\
            --out-dir . \\
            --coverage-files files.list \\
            ${pon_flag} \\
            --genome ${params.genome_build}

        mv *rds ${output}

        cp .command.out purecn_normaldb.log
        """
    }
    //
}
