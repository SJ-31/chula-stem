process DELLY {
    ext version: "1.3.1"
    conda { task.ext.conda }
    publishDir "$meta.out", mode: "copy"
    publishDir "$meta.log", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(tumor), path(normal)
    val(reference)
    val(exclude)
    val(module_number)
    //

    output:
    tuple val(meta), path(out)
    path("*.log")
    //

    script:
    out = "${module_number}-${meta.id}_Delly.bcf.gz" // TODO: you don't know what
    // the output of this is yet
    check = file("${meta.out}/${out}")
    if (check.exists()) {
        """
        cp ${check} .
        cp ${meta.log}/delly.log .
        """
    } else {
        """
        delly call \\
            -g ${reference} \\
            -x ${exclude} \\
            -o ${out} \\
            ${tumor} \\
            ${normal}

        cp .command.out delly.log
        """
    }
    //
}
