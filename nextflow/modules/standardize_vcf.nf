// Collect and standardize variant metrics from vcf files as a tsv
// Includes...
// - Variant Allele Fraction (INFO/AF): with `ChromosomeCounts`
// - Read depth (INFO/DP): with `Coverage`
// - Number of times alt allele is represented/allele count (INFO/AC): with `ChromosomeCounts`
// - Allele depth (FORMAT/AD): with `DepthPerAlleleBySample`
// - Median mapping quality (INFO/MMQ): with `MappingQuality`
// - Median base quality (INFO/MBQ): with `BaseQuality`
// - the total number of alleles in genotypes (INFO/AN): `ChromosomeCounts`
//
// Because some callers may use the above INFO or FORMAT tags for other purposes, they are first cleared before being (re)calculated
// Note: We rely on the callers for FORMAT/DP (SV callers don't produce this)

process STANDARDIZE_VCF {
    ext version: params.gatk_version

    publishDir "${meta.out}", mode:"copy", saveAs: params.saveFn
    publishDir "${meta.log}", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(vcf), path(normal_bam), path(tumor_bam)
    val(reference)
    val(module_number)

    output:
    tuple val(meta), path(output), emit: vcf
    path("*.log")
    //

    script:
    suffix = meta.suffix ? meta.suffix : "Standardized"
    output = "${module_number}-${meta.filename}-${suffix}.vcf.gz"
    check = file("${meta.out}/${output}")
    to_clear = [
        "INFO": ["DP", "MMQ", "MBQ", "AN", "AC", "AF"],
        "FORMAT": ["AD"],
    ]
    annotations = ["ChromosomeCounts", "Coverage",
                   "DepthPerAlleleBySample", "MappingQuality", "BaseQuality"]

    annotation_flag = annotations.collect({" -A ${it} "}).join(" ")
    clear_flag = to_clear.collect({ k,v -> v.collect({ "${k}/${it}" }).join(",") }).join(",")
    args = task.ext.args.join(" ")
    if (check.exists()) {
        """
        ln -sr ${check} .
        ln -sr ${meta.log}/${suffix}.log .
        """
    } else {
        """
        bcftools norm --fasta-ref ${reference} --atomize ${vcf} | \\
            bcftools annotate -x ${clear_flag} -O z > tmp.vcf.gz

        gatk IndexFeatureFile -I tmp.vcf.gz
        gatk VariantAnnotator -R ${reference} -V tmp.vcf.gz \\
            -I ${tumor_bam} -I ${normal_bam} \\
            ${annotation_flag} \\
            -O ${output} \\
            --add-output-vcf-command-line

        get_nextflow_log.bash ${suffix}.log
        """
    }
    //
}
