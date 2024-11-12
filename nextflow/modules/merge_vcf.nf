process MERGE_VCF {
    ext version: "1.21"
   

    publishDir "$meta.out", mode: "copy"
    publishDir "$meta.log", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(vcfs)
    val(module_number)
    //

    output:
    tuple val(meta), path(output), emit: vcf
    path("*.log")
    //

    script:
    output = "${module_number}-${meta.id}_merged.vcf.gz"
    check = file("${meta.out}/${output}")
    if (check.exists()) {
        """
        cp $check .
        cp ${meta.log}/merge_vcf.log .
        """
    } else {
        """
        bcftools merge *.vcf.gz -o $output
        cp .command.out merge_vcf.log
        """
    }
    //
}
