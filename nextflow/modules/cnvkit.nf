process CNVKIT {
    ext version: "0.9.11"

    publishDir "$meta.out", mode: "copy", saveAs: { x -> x ==~ /.*\.log/ ? null : x }
    publishDir "$meta.log", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(tumor)
    val(cnn_reference) // A copy number reference file (".cnn") created with the pooled normal samples
    val(omics_type) // hybrid, amplicon or wgs
    val(module_number)
    //

    output:
    tuple val(meta), path("${out}/")
    path(out)
    path("*.log")
    //

    script:
    out = "${module_number}-${meta.filename}-Cnvkit"
    check = file("${meta.out}/${out}")
    if (check.exists()) {
        """
        cp -r $check .
        ln -sr ${meta.log}/cnvkit.log .
        """
    } else {
        """
        cnvkit.py batch \\
                ${tumor} \\
            -r ${cnn_reference} \\
            -d ${out}

        get_nextflow_log.bash cnvkit.log
        """
    }
    //
}
