process GRIDSS {
    ext version: "2.13.2"

    publishDir "meta.out", mode:"copy", saveAs: { x -> x ==~ /.*\.log/ ? null : x }
    publishDir "meta.log", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(normal), path(tumor), path(indices, arity: "2")
    val(reference)
    val(blacklist) // list of genomic regions to exclude, gridss docs recommends using
    // ENCODE blacklists
    val(module_number)
    //

    output:
    tuple val(meta), path("${filtered}.gz")
    path("${output}.gz")
    path("*.log")
    //

    script:
    prefix = "${module_number}-${meta.filename}"
    output = "${prefix}-Gridss_all.vcf"
    filtered = "${prefix}-Gridss_confident.vcf"
    check1 = file("${meta.out}/${output}.gz")
    check2 = file("${meta.out}/${filtered}.gz")
    args = task.ext.args.join(" ")
    if (check1.exists() && check2.exists()) {
        """
        ln -sr ${check2} .
        ln -sr ${check1} .
        ln -sr ${meta.log}/gridss .
        """
    } else {
        """
        gridss \\
            --threads ${task.ext.threads} \\
            --reference ${reference} \\
            --output tmp.vcf \\
            ${normal} \\
            ${tumor}

        gridss_somatic_filter \\
            --input tmp.vcf \\
            --output "${filtered}" \\
            --fulloutput "${output}"

        bgzip ${filtered}
        bgzip ${output}

        get_nextflow_log.bash gridss.log
        """
    }
    //
    // TOOD: you can supply a pon with this
}
