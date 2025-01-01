process CROSS_REFERENCE {
    ext version: "1"
    // Cross-reference CNV or repeat data with those found in online databases
    // Takes output from ClassifyCNV and the get_msisensor.bash script and adds the
    // metadata columns "accession", "source", "p_overlap" and "ClinGen_report"

    publishDir "${meta.out}", mode:"copy", saveAs: params.saveFn
    publishDir "${meta.log}", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(cnv_or_msi)
    val(type) // Either "CNV" or "MSI"
    val(clingen)
    val(reference)
    val(module_number)
    //

    output:
    tuple val(meta), path(out)
    path("*.log")
    //

    script:
    out = Utl.getName(module_number, meta, "CR", "tsv")
    check = file("${meta.out}/${out}")
    if (check.exists()) {
        """
        ln -sr ${check} .
        ln -sr ${meta.log}/${type}_cross_reference .
        """
    } else {
        """
        cross_reference.r --input ${cnv_or_msi} \\
            --type ${type} \\
            --reference ${reference} \\
            --clingen ${clingen} \\
            --output ${out}

        cp .command.out ${type}_cross_reference.log
        """
    }
    //
}
