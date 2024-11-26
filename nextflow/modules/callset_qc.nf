process CALLSET_QC {
    // By default, filters a callset by minimum read depth in the tumor and maximum read depth in normal
    // as well as with the FILTER column, e.g. checking for 'PASS'
    ext version: "1.21"

    publishDir "${meta.out}", mode:"copy", saveAs: params.saveFn
    publishDir "${meta.log}", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(vcf)
    val(suffix)
    val(min_tumor_depth)
    val(max_normal_depth)
    val(accepted_filters) // List of values in the FILTER column that would permit a call
    val(module_number)
    //

    output:
    tuple val(meta), path(output)
    path("*.log")
    //

    script:
    suffix = suffix == "" ? "QC" : suffix
    output = "${module_number}-${meta.filename}-${suffix}.vcf.gz"
    check = file("${meta.out}/${output}")
    args = task.ext.args.join(" ")

    if (min_tumor_depth && max_normal_depth) {
        depth_flag = "-i 'FORMAT/DP[@normal.txt] <= ${max_normal_depth} && FORMAT/DP[@tumor.txt] >= ${min_tumor_depth}' " // SV callers don't have FORMAT/DP
    } else {
        depth_flag = ""
    }

    filter_flag = accepted_filters.collect({ "FILTER~'${it}'" }).join(" && ")

    if (check.exists()) {
        """
        ln -sr ${check} .
        ln -sr ${meta.log}/filter_qc.log .
        """
    } else {
        """
        echo ${meta.RGSM_normal} > normal.txt
        echo ${meta.RGSM_tumor} > tumor.txt

        bcftools filter -e 'FILTER="."' ${vcf} | \
            bcftools filter -i "${filter_flag}" | \
            bcftools filter ${depth_flag} -O z > ${output}

        get_nextflow_log.bash filter_qc.log
        """
    }
    //
}
