process CALLSET_QC_TSV {
    // Implements the same call filters as "CALLSET_QC", but has
    // additional formatting for VEP and options for merging results of different callers
    //
    ext version: "1"

    publishDir "${meta.out}", mode:"copy", saveAs: params.saveFn
    publishDir "${meta.log}", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(tsv)
    val(qc)
    // Meta must have a key "qc" whose value is a map specifying filters to apply
    // Currently supported are..
    // -- General
    //  - min_tumor_depth: the minimum number of ALT reads found in the tumor (FORMAT/DP)
    //  - max_normal_depth: the maximum number of ALT reads found in the normal (FORMAT/DP)
    //  - min_vaf: minimum variant allele frequency in the tumor (FORMAT/VAF)
    //  - accepted_filters: A string of FILTER flags that calls must have to be accepted,
    //      separated by ","
    //
    // -- Caller merging
    //  - min_callers: variants must be called by this number of different callers
    //  - VAF_adaptive: use VAF adaptive version of merging callers
    //
    // -- Resolving transcripts (flags)
    //  - canonical: keep only canonical transcripts whenever possible
    //  - informative: keep transcripts with the least number of NA VEP cols
    //  - impact: keep transcripts with the highest impact
    val(module_number)
    //

    output:
    tuple val(meta), path(output), emit: vcf
    path("*.log")
    //

    script:
    output = Utils.getName(module_number, meta, "QC", "tsv")
    check = file("${meta.out}/${output}")
    args = task.ext.args.join(" ")

    flags = ["informative", "canonical", "impact", "vaf_adaptive"]
    qc_copy = qc ? qc.clone() : meta.qc.clone()
    accepted_filters = qc_copy.remove("accepted_filters")
    if (accepted_filters) {
        joined = accepted_filters.join(",")
        filter_flag = "--accepted_filters ${joined}"
    } else {
        filter_flag = ""
    }
    all = qc_copy.collect({ k, v ->
        if (k in flags) {
            v ? "--${k}" : ""
        } else {
            v ? "--${k} ${v}" : ""
        }
    }).join(" ")

    if (check.exists()) {
        """
        ln -sr ${check} .
        ln -sr ${meta.log}/filter_qc.log .
        """
    } else {
        """
        callset_qc -i ${tsv} \\
            --tool_source_tag ${params.source_description} \\
            -o ${output} \\
            ${args} \\
            ${filter_flag} \\
            ${all}

        get_nextflow_log.bash qc_tsv.log
        """
    }
    //
}
