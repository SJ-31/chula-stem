process CNVKIT {
    ext version: "0.9.11"

    publishDir "$meta.out", mode: "copy", saveAs: params.saveFn
    publishDir "$meta.log", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(tumor), path(snps), path(purity), path(ploidy)
    tuple val(cnn_reference)
    // A copy number reference file (".cnn") created with the pooled normal samples
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
    purity_val = file(purity).read // TODO: don't know if this works correctly
    ploidy_val = file(ploidy).read
    if (check.exists()) {
        """
        cp -r ${check} .
        ln -sr ${meta.log}/cnvkit.log .
        """
    } else {
        """
        cnvkit.py batch \\
                ${tumor} \\
            -r ${cnn_reference} \\
            -d ${out}

        cnvkit.py call \\
            ${out}/${out}.call.cns \\
            --vcf ${snps} \\
            --purity ${purity_val} \\
            --ploidy ${ploidy_val} \\
            --method clonal \\
            --output ${out}.call.clonal.cns

        get_nextflow_log.bash cnvkit.log
        """
    }
    // Will perform a second call to adjust original with information of tumor purity, ploidy
    // and allele frequencies from known variants
    // The or
}
