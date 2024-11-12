process SNPEFF {
    ext version: "5.2e"
    

    publishDir "$meta.out", mode: "copy"
    publishDir "$meta.log", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(vcf)
    val(reference)
    val(sample_definition)
    val(module_number)
    //

    output:
    path("${output}.gz")
    path("*.html")
    path("*.log")
    //

    script:
    output = "${module_number}-${meta.id}_snpEff.vcf"
    report = "${module_number}-${meta.id}_snpEff_summary.html"
    check = file("$meta.out/${output}.gz")
    if (check.exists()) {
        """
        cp $check .
        cp "${meta.out}/${report}" .
        cp "${meta.log}/snpEff.log" .
        """
    } else {
        """
        $params.snpEff -cancer \\
            -canon \\
            -nodownload  \\
            -cancerSamples $sample_definition \\
            -v \\
            $reference \\
            $vcf > $output
        gzip $output

        cp snpEff_summary.html $report
        cp .command.out snpEff.log
        """
    }
    //
}
