process CALLSET_QC {
    // By default, filters a callset by minimum read depth in the tumor and maximum read depth in normal
    // as well as with the FILTER column, e.g. checking for 'PASS'
    ext version: "1.21"

    publishDir "${meta.out}", mode:"copy", saveAs: params.saveFn
    publishDir "${meta.log}", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(vcf)
    val(qc) // A map specifying filters to apply
    // Currently supported are..
    // - min_tumor_depth: the minimum number of ALT reads found in the tumor (FORMAT/DP)
    // - max_normal_depth: the maximum number of ALT reads found in the normal (FORMAT/DP)
    // - min_VAF: minimum variant allele frequency in the tumor (FORMAT/VAF)
    // - accepted_filters: A list of FILTER flags that calls must have to be accepted
    val(panel_of_normals)
    val(module_number)
    //

    output:
    tuple val(meta), path(output), emit: vcf
    path("*.log")
    //

    script:
    output = Utl.getName(module_number, meta, "QC", "vcf.gz")
    check = file("${meta.out}/${output}")
    args = task.ext.args.join(" ")

    ndepth = qc.max_normal_depth ? "FORMAT/AD[@normal.txt:1-] <= ${qc.max_normal_depth}" : ""
    tdepth = qc.min_tumor_depth ? "FORMAT/AD[@tumor.txt:1-] >= ${qc.min_tumor_depth}" : ""
    vaf = qc.min_VAF ? "FORMAT/AF[@tumor.txt] > ${qc.min_VAF}" : ""
    // Note: VAF is calculated with GATK, which calls it FORMAT/AF
    filter = qc.accepted_filters.collect({ "FILTER~\"${it}\"" }).join(" && ")

    filter_list = [tdepth, filter, vaf]
    if (!params.tumor_only) {
        filter_list << ndepth
    }
    all = filter_list.findAll({ it != "" && it != null })
        .collect({ "bcftools filter -i '${it}'" }).join(" | ")
    all = all != "" && all != null ? "${all} -O z > ${output}" : "bcftools view -O z > ${output}"

    input = panel_of_normals ? "pon_filtered.vcf.gz" : vcf
    if (panel_of_normals) {
        pon_filter_command = """
        bcftools view -T ^${panel_of_normals} ${vcf} -O z -W -o pon_filtered.vcf.gz 
        """
    } else {
        pon_filter_command = ""
    }

    if (check.exists()) {
        """
        ln -sr ${check} .
        ln -sr ${meta.log}/filter_qc.log .
        """
    } else {
        """
        ${pon_filter_command}

        echo ${meta.RGSM_normal} > normal.txt
        echo ${meta.RGSM_tumor} > tumor.txt

        bcftools filter -e 'FILTER="."' ${input} | ${all}

        get_nextflow_log.bash filter_qc.log
        """
    }
    //
}
