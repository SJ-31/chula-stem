process FACETS {
    ext version: "2"

    publishDir "$meta.out", mode: "copy", saveAs: params.saveFn
    publishDir "$meta.log", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(counts_file), val(window_size) // Output from facets_pileup
    val(module_number)
    //

    output:
    tuple val(meta), path(prefix), emit: cnv
    tuple val(with_caller), path("${prefix}/${prefix}_hisens.rds"), emit: rds
    tuple val(meta), path("${prefix}/purity.txt"), path("${prefix}/ploidy.txt"), emit: purity_ploidy
    path(tsv)
    path("*.log")
    //

    script:
    prefix = Utl.getName(module_number, meta, "Facets")
    check = file("${meta.out}/${prefix}")
    tsv = "${prefix}/${prefix}_hisens.tsv"
    check2 = file("${meta.out}/${tsv}")
    with_caller = meta + ["caller": "facets"]
    window_size = window_size ? window_size : 250
    args = task.ext.args.join(" ")
    if (check.exists() && check2.exists()) {
        """
        cp -r $check .
        ln -sr ${check2} .
        ln -sr ${meta.log}/facets.log .
        """
    } else {
        """
        run-facets-wrapper.R \\
            --counts-file ${counts_file} \\
            --sample-id ${prefix} \\
            --genome ${task.ext.genome} \\
            --snp-window-size ${window_size} \\
            ${args}

        cut -f 2 ${prefix}/${prefix}.txt | tail -n 1 > ${prefix}/purity.txt
        cut -f 3 ${prefix}/${prefix}.txt | tail -n 1 > ${prefix}/ploidy.txt

        get_facets.bash "${prefix}/${prefix}_hisens.rds" "${tsv}"

        get_nextflow_log.bash facets.log
        """
    }
    //
}
