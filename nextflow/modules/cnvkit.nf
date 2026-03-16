process CNVKIT {
    ext version: "0.9.11"
    errorStrategy "ignore" // [2025-06-24 Tue] Temporary incompatibility with cnvkit
    // and the pdac reference file

    // BUG: cnvkit doesn't like the GT format header and having multiple samples
    // if the data are paired, should get only the tumor sample and clean it up
    
    publishDir "$meta.out", mode: "copy", saveAs: params.saveFn
    publishDir "$meta.log", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(tumor), path(snps), val(purity), val(ploidy)
    val(cnn_reference)
    // A copy number reference file (".cnn") created with the pooled normal samples
    val(omics_type) // hybrid, amplicon or wgs
    val(module_number)
    //

    output:
    tuple val(with_caller), path("${out}/${out}.call.cns"), emit: cns
    tuple val(with_caller), path("${out}/${out}.cnr"), emit: cnr
    path(out), emit: all
    path("*.log")
    //

    script:
    out = Utl.getName(module_number, meta, "Cnvkit")
    check = file("${meta.out}/${out}")
    with_caller = meta + ["caller": "cnvkit"]
    ploidy_val = ploidy ? ploidy : params.defaults.ploidy
    purity_val = purity ? purity : params.defaults.purity
    name = tumor.baseName

    clonal = !ploidy_val && !purity_val ? "false" : "true"
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

        cd "${out}"
        for file in ${name}*; do
            prefix=\$(echo "\${file}" | sed 's/${name}//')
            mv \${file} "${out}\${prefix}"
        done
        cd ..

        if [[ "${clonal}" == "true" ]]; then
            cnvkit.py call \\
                "${out}/${out}.call.cns" \\
                --vcf "${snps}" \\
                --purity "${purity_val}" \\
                --ploidy "${ploidy_val}" \\
                --method clonal \\
                --output "${out}.call.clonal.cns"
            mv "${out}.call.clonal.cns" "${out}"
        fi

        get_nextflow_log.bash cnvkit.log
        """
    }
    // Will perform a second call to adjust original with information of tumor purity, ploidy
    // and allele frequencies from known variants
}
