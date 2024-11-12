process MUTECT2 {
    ext version: "4.6.1.0"

    publishDir "$meta.out", mode: "copy"
    publishDir "$meta.log", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(tumor), path(normal)
    val(reference)
    val(module_number)
    //

    output:
    tuple val(meta), path(out), emit: variants
    path("*.log")
    //

    script:
    out = "${module_number}-${meta.id}_Mutect2.vcf.gz"
    check = file("${meta.out}/${out}")
    if (check.exists()) {
        """
        cp $check .
        cp ${meta.log}/mutect2.log .
        """
    } else {
        """
        gatk Mutect2 \\
            -R $reference \\
            -I $tumor \\
            -I $normal \\
            -normal $meta.RGSM_normal \\
            --output temp.vcf.gz

        vcf_info_add_tag -n SOURCE \\
            -d "Tool producing call" \\
            -b '.' \\
            -t String \\
            -a mutect2 \\
            -i temp.vcf.gz \\
            -o $out

        cp .command.out mutect2.log
        """
    }
    //
}
// Look into using --panel-of-normals flag \\
