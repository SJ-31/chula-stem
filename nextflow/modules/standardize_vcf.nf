// Collect and standardize variant metrics from vcf files as a tsv
// Includes...
// - Variant Allele Fraction (FORMAT/AF): with `ChromosomeCounts`
// - Read depth (INFO/DP): with `Coverage`
// - Number of times alt allele is represented/allele count (INFO/AC): with `ChromosomeCounts`
// - Allele depth (AD): with `DepthPerAlleleBySample`
// - Median mapping quality (INFO/MMQ): with `MappingQuality`
// - Median base quality (INFO/MBQ): with `BaseQuality`
// - the total number of alleles in genotypes (AN): `ChromosomeCounts`
//
// Because some callers may use the above INFO or FORMAT tags for other purposes, they are first cleared before being (re)calculated

process STANDARDIZE_VCF {
    ext version: params.gatk_version

    publishDir "${meta.out}", mode:"copy", saveAs: params.saveFn
    publishDir "${meta.log}", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(vcf), path(normal_bam), path(tumor_bam)
    val(reference)
    val(target_intervals)
    val(module_number)

    output:
    tuple val(meta), path(output), emit: vcf
    path("*.log")
    //

    script:
    output = "${module_number}-${meta.filename}-standardized.vcf.gz"
    check = file("${meta.out}/${output}")
    to_clear = [
        "INFO": ["DP", "MMQ", "MBQ", "AN", "AC"],
        "FORMAT": ["AF", "AN", "AD", "DP"],
    ]
    annotations = ["ChromosomeCounts", "Coverage",
                   "DepthPerAlleleBySample", "MappingQuality", "BaseQuality"]

    annotation_flag = annotations.collect({" -A ${it} "}).join(" ")
    clear_flag = to_clear.collect({ k,v -> v.collect({ "${k}/${it}" }).join(",") }).join(",")
    target_flag = target_intervals != "" ? "--intervals ${target_intervals}" : ""
    args = task.ext.args.join(" ")
    if (check.exists()) {
        """
        cp -r ${check} .
        ln -sr ${meta.log}/bcftools_stats.log .
        """
    } else {
        """
        bcftools annotate -x ${clear_flag} ${vcf} -O z > tmp.vcf.gz
        gatk IndexFeatureFile -I tmp.vcf.gz
        gatk VariantAnnotator ${target_flag} -R ${reference} -V tmp.vcf.gz \\
            -I ${tumor_bam} -I ${normal_bam} \\
            ${annotation_flag} \\
            -O ${output}

        get_nextflow_log.bash standardize_vcf.log
        """
    }
    //
}
