process MUTECT2 {
    publishDir "$meta.out", mode: "copy"
    publishDir "$meta.log", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(tumor), path(normal)
    val(reference)
    val(module_number)
    //

    output:
    tuple val(meta), path(out)
    path("*.log")
    //

    script:
    out = "${module_number}-${meta.baseName}_Mutect2.vcf.gz"
    def check = file("${meta.out}/${out}")
    if (check.exists()) {
        """
        cp $check.name .
        cp ${meta.log}/mutect2.log .
        """
    } else {
        """
        gatk Mutect2 \
            -R $reference \
            -I $tumor \
            -I $normal \
            -normal $normal.baseName \
            --panel-of-normals ??? \
            --output temp.vcf.gz

        bcftools annotate temp.vcf.gz TODO
        cp .command.out mutect2.log
        """
    }
    //
}
