process MUTECT2 {
    ext version: "4.6.1.0"

    publishDir "$meta.out", mode: "copy", saveAs: { x -> x ==~ /.*\.log/ ? null : x }
    publishDir "$meta.log", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(normal), path(tumor), path(indices, arity: "2")
    val(reference)
    val(module_number)
    //

    output:
    tuple val(meta), path(out), emit: variants
    path(stats)
    path("*.log")
    //

    script:
    out = "${module_number}-${meta.filename}-Mutect2.vcf.gz"
    stats = "${module_number}-${meta.filename}-Mutect2_stats.txt"
    uncompressed = out.replace(".gz", "")
    check = file("${meta.out}/${out}")
    if (check.exists()) {
        """
        ln -sr $check .
        ln -sr ${meta.log}/mutect2.log .
        ln -sr "${meta.out}/${stats}" .
        """
    } else {
        """
        gatk Mutect2 \\
            -R $reference \\
            -I $tumor \\
            -I $normal \\
            -normal $meta.RGSM_normal \\
            --output temp.vcf.gz

        mv temp.vcf.gz.stats $stats

        vcf_info_add_tag -n SOURCE \\
            -d "$params.source_description" \\
            -b '.' \\
            -t String \\
            -a mutect2 \\
            -i temp.vcf.gz \\
            -o $uncompressed

        bgzip $uncompressed
        get_nextflow_log.bash mutect2.log
        """
    }
    //
}
// Look into using --panel-of-normals flag \\
