process FACETS_PILEUP {
    ext version: "2"

    publishDir "$meta.out", mode: "copy", saveAs: params.saveFn
    publishDir "$meta.log", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(normal), path(tumor), path(indices, arity: "2")
    val(snps) // "should contain polymorphic SNPs, so that FACETS can infer changes in allelic configuration at genomic loci from changes in allele ratios"
    // Docs recommend dbSNP
    val(module_number)
    //

    output:
    tuple val(meta), path(output), emit: pileup
    path("*.log")
    //

    script:
    prefix = params.getName(module_number, meta, "Facets_pileup")
    output = "${prefix}.snp_pileup.gz"
    check = file("${meta.out}/${output}")
    if (check.exists()) {
        """
        ln -sr $check .
        ln -sr ${meta.log}/facets_pileup.log .
        """
    } else {
        """
        snp-pileup-wrapper.R \\
            --vcf-file ${snps} \\
            --normal-bam ${normal} \\
            --tumor-bam ${tumor} \\
            --output-prefix ${prefix}

        get_nextflow_log.bash facets_pileup.log
        """
    }
    //
}
