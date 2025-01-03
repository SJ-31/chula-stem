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
    errorStrategy "ignore"
    // BUG: <2025-01-03 Fri> gatk has an Invalid File pointer error
    // when trying to read the recal bam file it produced...

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
    output = Utl.getName(module_number, meta, "Standardize", "vcf.gz")
    check = file("${meta.out}/${output}")
    to_clear = [
        "INFO": ["DP", "MMQ", "MBQ", "AN", "AC", "AF"],
        "FORMAT": ["AD"],
    ]
    annotations = ["ChromosomeCounts", "Coverage", "AlleleFraction",
                   "DepthPerAlleleBySample", "MappingQuality", "BaseQuality"]

    annotation_flag = annotations.collect({" -A ${it} "}).join(" ")
    clear_flag = to_clear.collect({ k,v -> v.collect({ "${k}/${it}" }).join(",") }).join(",")
    read_flag = !params.tumor_only ? " -I ${tumor_bam} -I ${normal_bam} " : " -I ${tumor_bam} "
    args = task.ext.args.join(" ")
    if (check.exists()) {
        """
        ln -sr ${check} .
        ln -sr ${meta.log}/${suffix}.log .
        """
    } else {
        """
        standardize_vcf_clean.bash "${clear_flag}" ${vcf}

        gatk IndexFeatureFile -I tmp.vcf.gz
        gatk VariantAnnotator ${args} \\
            -R ${reference} -V tmp.vcf.gz \\
            ${read_flag} \\
            ${annotation_flag} \\
            -O ${output} \\
            --add-output-vcf-command-line

        get_nextflow_log.bash ${suffix}.log
        """
    }
    //
}
