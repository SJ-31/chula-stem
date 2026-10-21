process MSISENSORPRO {
    ext version: "1.3.0"

    input:
    tuple val(meta), path(normal), path(tumor), path(indices, arity: "2")
    val reference // A homopolymers and microsatellites tsv file, generated
    // with msisensor-pro scan -d <reference genome> -o <tsv>
    val omics_type // Either exome or wgs
    val gff // Path to genome gff file to identify genes in repetitive regions

    output:
    tuple val(meta), path("${prefix}/${unstable}.tsv"), emit: tsv
    tuple val(meta), path("${prefix}/${prefix}_summary.tsv"), emit: summary
    tuple val(meta), path(prefix), emit: all
    tuple val(meta), path("*.log"), emit: log

    script:
    //
    prefix = Utl.getName(meta, "Msisensor") // Is the prefix, but also the summary file
    all = "${prefix}_all"
    unstable = "${prefix}_unstable" // Sites in "all" with statistically significant p-values
    distribution_file = "${prefix}_dis"

    command = !params.tumor_only ? "msi" : "pro"
    n_flag = !params.tumor_only ? " -n ${normal} " : ""
    cov_threshold = omics_type == "exome" ? 20 : 15

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
        -s "${prefix}"_summary.tsv \\
        -g ${gff}

    mkdir ${prefix}
    mv *tsv ${prefix}
    get_nextflow_log.bash msisensor.log
    """
}

process MSISENSORPRO_COLLECT {
    input:
    tuple val(meta), path(summaries)

    output:
    path output

    script:
    output = Utl.getName(meta, "Msisensor_summaries", "tsv")
    """
    collect_msi.py -i . \
        -s Msisensor_summary.tsv \
        -o ${output}
    """
}
