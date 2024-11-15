process CNVKIT {
    ext version: "0.9.11"

    publishDir "$meta.out", mode: "copy", saveAs: { x -> x ==~ /.*\.log/ ? null : x }
    publishDir "$meta.log", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(normal), path(tumor), path(indices, arity: "2")
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
        ln -sr $check .
        ln -sr ${meta.log}/cnvkit.log .
        """
    } else {
        """
        cnvkit.py batch \\
            --normal ${normal} \\
            -r ${reference} \\
            -d ${out} \\
            ${tumor}

        get_nextflow_log.bash msisensor.log
        """
    }
    //
}
