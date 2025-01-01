process DELLY_COV {
    ext version: "1.3.1"

    publishDir "$meta.out", mode: "copy", saveAs: params.saveFn
    publishDir "$meta.log", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(bam), path(indices, arity: "1..*")
    val(reference)
    val(mappability)
    val(module_number)
    //

    output:
    tuple val(meta), path(output), emit: cov
    path("*.log")
    //

    script:
    output = Utl.getName(module_number, meta, "Delly", "cov.gz")
    check = file("${meta.out}/${output}")
    args = task.ext.args.join(" ")
    if (check.exists()) {
        """
        ln -sr ${check} .
        ln -sr ${meta.log}/delly_cov.log .
        """
    } else {
        """
        delly cnv -a \\
            --genome ${reference} \\
            --mappability ${mappability} \\
            --covfile ${output} \\
            --outfile tmp.bcf \\
            ${bam}

        get_nextflow_log.bash delly_cov.log
        """
    }
    //
}
