process MSISENSORPRO {
    ext version: "1.3.0"
    conda { task.ext.conda }

    publishDir "$meta.out", mode: "copy"
    publishDir "$meta.log", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(tumor), path(normal)
    val(reference) // A homopolymers and microsatellites tsv file, generated
    // with msisensor-pro scan -d <reference genome> -o <tsv>
    val(module_number)
    //

    output:
    tuple val(meta), path(out)
    path("*.log")
    //

    script:
    out = "${module_number}-${meta.id}_Msisensor.vcf.gz" // TODO: you don't know what
    // the output of this is yet
    check = file("${meta.out}/${out}")
    if (check.exists()) {
        """
        cp ${check} .
        cp ${meta.log}/msisensor.log .
        """
    } else {
        """
        msisensor-pro msi \\
            -n ${normal} \\
            -t ${tumor} \\
            -d ${reference} \\
            -o ${out}

        cp .command.out msisensor.log
        """
    }
    //
}
