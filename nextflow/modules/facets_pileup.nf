process FACETS_PILEUP {
    ext version: "2"

    input:
    tuple val(meta), path(normal), path(tumor), path(indices, arity: "2")
    val snps // "should contain polymorphic SNPs, so that FACETS can infer changes in allelic configuration at genomic loci from changes in allele ratios"

    output:
    // Docs recommend dbSNP

    //
    tuple val(meta), path(output), emit: pileup
    tuple val(meta), path("*.log"), emit: logfile

    script:
    //
    prefix = Utl.getName(meta, "Facets")
    output = "${prefix}.snp_pileup.gz"
    """
    snp-pileup-wrapper.R \\
        --vcf-file ${snps} \\
        --normal-bam ${normal} \\
        --tumor-bam ${tumor} \\
        --output-prefix ${prefix}

    get_nextflow_log.bash facets_pileup.log
    """
}
