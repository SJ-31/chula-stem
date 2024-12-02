process CNVKIT {
    ext version: "0.9.11"

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
    tuple val(with_caller), path("${out}/${out}.call.cns"), emit: cnv
    path(out)
    path("*.log")
    //

    shell:
    out = Utils.getName(module_number, meta, "Cnvkit")
    check = file("${meta.out}/${out}")
    with_caller = meta + ["caller": "cnvkit"]
    ploidy_val = ploidy ? ploidy : params.defaults.ploidy
    purity_val = purity ? purity : params.defaults.purity

    clonal = !ploidy_val && !purity_val ? "false" : "true"
    if (check.exists()) {
        '''
        cp -r !{check} .
        ln -sr !{meta.log}/cnvkit.log .
        '''
    } else {
        template 'cnvkit.bash'
    }
    // Will perform a second call to adjust original with information of tumor purity, ploidy
    // and allele frequencies from known variants
}
