process CALLSET_QC_TSV {
    // Implements the same call filters as "CALLSET_QC", but has
    // additional formatting for VEP and options for merging results of different callers
    //
    ext version: "1"

    publishDir "${meta.out}", mode:"copy", saveAs: params.saveFn
    publishDir "${meta.log}", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(tsv), path(ignore_regions)
    // Will remove variants in locations that were called by the repetitive
    //      element tools, as those are likely to be more accurate in such cases
    //      Can be disabled by just sending an empty file
    val(qc)
    // Meta can also have a key "qc" whose value is a map specifying filters to apply
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
    tuple val(meta), path(output), emit: tsv
    path("*.log")
    //

    script:
    output = Utl.getName(module_number, meta, "QC", "tsv")
    check = file("${meta.out}/${output}")
    args = task.ext.args.join(" ")
    outlog = "${output}_qc_tsv.log"

    ignore_flag = !ignore_regions.empty() ? " --ignore_regions ${ignore_regions} " : ""

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
        ln -sr ${meta.log}/${outlog} .
        """
    } else {
        """
        callset_qc -i ${tsv} \\
            ${ignore_flag} \\
            --tool_source_tag "${params.source_name}" \\
            -o ${output} \\
            ${args} \\
            ${filter_flag} \\
            ${all}

        get_nextflow_log.bash ${outlog}
        echo "QC flags\n${all}" >> ${outlog}
        """
    }
    //
}
