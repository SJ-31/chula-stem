process FACETS {
    ext version: "2"

    publishDir "$meta.out", mode: "copy", saveAs: { x -> x ==~ /.*\.log/ ? null : x }
    publishDir "$meta.log", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(counts_file) // Output from facets_pileup
    val(module_number)
    //

    output:
    tuple path(meta), path(prefix)
    path(prefix)
    path("*.log")
    //

    script:
    prefix = "${module_number}-${meta.filename}-Facets"
    check = file("${meta.out}/${prefix}")
    args = task.ext.args.join(" ")
    if (check.exists()) {
        """
        cp -r $check .
        ln -sr ${meta.log}/facets.log .
        """
    } else {
        """
        run-facets-wrapper.R \\
            --counts-file ${counts_file} \\
            --sample-id ${prefix} \\
            ${args}

        get_nextflow_log.bash facets.log
        """
    }
    //
}
