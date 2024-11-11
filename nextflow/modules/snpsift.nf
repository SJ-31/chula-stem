process SNPSIFT {
    ext version: "5.2e"
    conda { task.ext.conda }

    publishDir "$meta.out", mode: "copy"
    publishDir "$meta.log", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(vcf)
    val(known_variants) // List of known somatic variants e.g.
    val(sample_definition)
    val(module_number)
    //

    output:
    path("${output}.gz")
    path("*.html")
    path("*.log")
    //

    script:
    def annotate_expr = known_variants.join(" ")
    output = "${module_number}-${meta.id}_SnpSift.vcf"
    check = file("$meta.out/${output}.gz")
    if (check.exists()) {
        """
        cp $check .
        cp "${meta.out}/${report}" .
        cp "${meta.log}/SnpSift.log" .
        """
    } else {
        """
        bcftools index $vcf

        $params.SnpSift filter $params.SnpfSift_filter $vcf | \\
            $params.SnpSift annotate $annotate_expr > \\
            $output

        gzip $output

        cp .command.out SnpSift.log
        """
    }
    //
}
