process CNVKIT {
    conda params.cnvkit_env
    // Version 0.9.11
    publishDir "$meta.out", mode: "copy"
    publishDir "$meta.log", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(tumor), path(normal)
    val(reference) // A copy number reference file (".cnn") created using
    // cnvkit.py reference -f <genome_fasta> -o <output>
    val(module_number)
    //

    output:
    path("*.log")
    //

    script:
    out = "${module_number}-${meta.id}_Cnvkit"
    def check = file("${meta.out}/${out}")
    if (check.exists()) {
        """
        cp ${check.name} .
        cp ${meta.log}/cnvkit.log .
        """
    } else {
        """
        cnvkit.py batch \\
            --normal ${normal} \\
            -r ${reference} \\
            -d ${out} \\
            ${tumor}

        cp .command.out msisensor.log
        """
    }
    //
}
