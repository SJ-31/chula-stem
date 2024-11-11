process CNVKIT {
    ext version: "0.9.11"
    conda { task.ext.conda }
    publishDir "$meta.out", mode: "copy"
    publishDir "$meta.log", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(tumor), path(normal)
    val(reference) // A copy number reference file (".cnn") created using
    // cnvkit.py reference -f <genome_fasta> -o <output>
    //  Note: if not using WGS, you should prepare a BED file listing genomic coordinates
    //      of captured regions during sample preparation. Might be able to find from vendor
    val(module_number)
    //

    output:
    path("*.log")
    //

    script:
    out = "${module_number}-${meta.id}_Cnvkit"
    check = file("${meta.out}/${out}")
    if (check.exists()) {
        """
        cp $check .
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
