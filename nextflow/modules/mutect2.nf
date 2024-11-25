process MUTECT2 {
    ext version: params.gatk_version

    publishDir "$meta.out", mode: "copy", saveAs: { x -> x ==~ /.*\.log/ ? null : x }
    publishDir "$meta.log", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(normal), path(tumor), path(indices, arity: "2")
    val(reference)
    val(target_intervals) // For exome data, path to target intervals file
    val(module_number)
    //

    output:
    tuple val(meta), path(out), emit: variants
    tuple val(meta), path(raw), emit: raw
    tuple val(meta), path(stats), emit: stats
    path("*.log")
    //

    script:
    out = "${module_number}-${meta.filename}-Mutect2.vcf.gz"
    stats = "${out}.stats"
    raw = "${module_number}-${meta.filename}-Mutect2_raw.tar.gz"
    target_flag = target_intervals == "" ? " --intervals ${target_intervals} " : ""
    uncompressed = out.replace(".gz", "")
    check = file("${meta.out}/${out}")
    if (check.exists()) {
        """
        ln -sr $check .
        ln -sr ${meta.log}/mutect2.log .
        ln -sr "${meta.out}/${stats}" .
        ln -sr "${meta.out}/${raw}" .
        """
    } else {
        """
        gatk Mutect2 \\
            -R $reference \\
            -I $tumor \\
            -I $normal \\
            -normal $meta.RGSM_normal \\
            ${target_flag} \\
            --f1r2-tar-gz $raw \\
            --output tmp.vcf.gz

        mv tmp.vcf.gz.stats $stats

        bcftools view -s "${meta.RGSM_normal},${meta.RGSM_tumor}" -O z tmp.vcf.gz > tmp2.vcf.gz

        vcf_info_add_tag -n SOURCE \\
            -d "$params.source_description" \\
            -b '.' \\
            -t String \\
            -a mutect2 \\
            -i tmp2.vcf.gz \\
            -o $uncompressed

        bgzip $uncompressed
        get_nextflow_log.bash mutect2.log
        """
    }
    //
}
// Look into using --panel-of-normals flag \\
// https://gatk.broadinstitute.org/hc/en-us/articles/360035531132--How-to-Call-somatic-mutations-using-GATK4-Mutect2 describes how to create a pon from your own data
