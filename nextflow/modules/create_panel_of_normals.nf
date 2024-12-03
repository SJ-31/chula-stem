process CREATE_PANEL_OF_NORMALS {
    ext version: params.gatk_version

    publishDir "${meta.out}", mode:"copy", saveAs: params.saveFn
    publishDir "${meta.log}", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(genomics_db)
    val(reference)
    val(module_number)
    //

    output:
    tuple val(meta), path(output), emit: pon
    path("*.log")
    //

    script:
    output = Utils.getName(module_number, meta, "PON", "vcf.gz")
    check = file("${meta.out}/${output}")
    args = task.ext.args.join(" ")
    if (check.exists()) {
        """
        ln -sr ${check} .
        ln -sr ${meta.log}/create_panel_of_normals.log .
        """
    } else {
        """
        gatk CreateSomaticPanelOfNormals ${args} \\
            -R ${reference} -V ${genomics_db} \\
            -O ${output}

        get_nextflow_log.bash create_panel_of_normals.log
        """
    }
    //
}
