process MSISENSORPRO {
    ext version: "1.3.0"

    publishDir "$meta.out", mode: "copy", saveAs: params.saveFn
    publishDir "$meta.log", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(normal), path(tumor), path(indices, arity: "2")
    val(reference) // A homopolymers and microsatellites tsv file, generated
    // with msisensor-pro scan -d <reference genome> -o <tsv>
    val(omics_type) // Either exome or wgs
    val(gff) // Path to genome gff file to identify genes in repetitive regions
    val(module_number)
    //

    output:
    tuple val(meta), path("*"), emit: tsvs
    path("*.log")
    //

    script:
    prefix = Utils.getName(module_number, meta, "Msisensor") // Is the prefix, but also the summary file
    all = "${prefix}_all"
    unstable = "${prefix}_unstable" // Sites in "all" with statistically significant p-values
    distribution_file = "${prefix}_dis"

    command = !params.tumor_only ? "msi" : "pro"
    n_flag = !params.tumor_only ? " -n ${normal} " : ""
    cov_threshold = omics_type == "exome" ? 20 : 15

    check_summary = file("${meta.out}/${prefix}_summary.tsv")
    check_all = file("${meta.out}/${all}.tsv")
    check_unstable = file("${meta.out}/${unstable}.tsv")
    check_dis = file("${meta.out}/${distribution_file}.txt")

    if (check_all.exists() && check_summary.exists() && check_dis.exists() && check_unstable.exists()) {
        """
        ln -sr $check_all .
        ln -sr $check_unstable .
        ln -sr $check_dis .
        ln -sr $check_summary .
        ln -sr ${meta.log}/msisensor.log .
        """
    } else {
        """
        msisensor-pro ${command} \\
            ${n_flag} \\
            -t ${tumor} \\
            -d ${reference} \\
            -o ${prefix} \\
            -c ${cov_threshold}

        mv "${prefix}" "${prefix}"_summary.tsv
        mv "${all}" "${all}".tsv
        mv "${distribution_file}" "${distribution_file}".txt

        get_msisensor.bash -i ${unstable} \\
            -o "${unstable}".tsv \\
            g ${gff}

        get_nextflow_log.bash msisensor.log
        """
    }
    //
}
