process CROSS_REFERENCE {
    ext version: "1"

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
    path("*.log")
    //

    script:
    out = Utils.getName(module_number, meta, "CR", "tsv")
    check = file("${meta.out}/${out}")
    if (check.exists()) {
        """
        ln -sr ${check} .
        ln -sr ${meta.log}/ .
        """
    } else {
        """
        cross_reference.r --input ${cnv_or_msi} \\
            --type ${type} \\
            --reference ${reference} \\
            --clingen ${clingen}

        cp .command.out .log
        """
    }
    //
}
