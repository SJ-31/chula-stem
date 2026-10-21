process SNPSIFT {
    ext version: "5.2e"
  
    
    

    input:
    tuple val(meta), path(vcf)
    val(known_variants) // List of known somatic variants e.g.
    val(sample_definition)
    
    //

    output:
    path("${output}.gz")
    path("*.html")
    tuple val(meta), path("*.log"), emit: log
    //

    script:
    def annotate_expr = known_variants.join(" ")
    output = Utl.getName(meta, "SnpSift", "vcf")
    check = file("$meta.out/${output}.gz")
    if (check.exists()) {
        """
        ln -sr $check .
        ln -sr "${meta.out}/${report}" .
        ln -sr "${meta.log}/SnpSift.log" .
        """
    } else {
        """
        bcftools index $vcf

        $params.SnpSift filter ${task.ext.args} $vcf | \\
            $params.SnpSift annotate $annotate_expr | \\
            $params.SnpSift varType > \\
            $output

        bgzip $output

        get_nextflow_log.bash SnpSift.log
        """
    }
    //
}
