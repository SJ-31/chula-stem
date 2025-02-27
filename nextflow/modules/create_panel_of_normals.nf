process CREATE_PANEL_OF_NORMALS {
    ext version: "1"

    publishDir "${meta.out}", mode:"copy", saveAs: params.saveFn
    publishDir "${meta.log}", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), val(vcfs)
    val(minimum)
    val(other_vcfs)
    val(module_number)
    //

    output:
    tuple val(meta), path(output), emit: pon // Empty pon with no sample information
    tuple val(meta), path(output_ws), emit: pon_ws
    path("*.log")
    //

    script:
    output = Utl.getName(module_number, meta, "PON", "vcf.gz")
    output_ws = Utl.getName(module_number, meta, "PON_ws", "vcf.gz")
    to_sample_spec = vcfs.toList().join("\n")
    to_other_spec = other_vcfs.join("\n")
    check = file("${meta.out}/${output}")
    if (check.exists()) {
        """
        ln -sr ${check} .
        ln -sr ${meta.log}/create_panel_of_normals.log .
        """
    } else {
        """
        echo -e -n "${to_sample_spec}" > samples.txt
        echo -e -n "${to_other_spec}" > other.txt

        create_panel_of_normals.bash -l samples.txt \\
            -m ${minimum} \\
            -v other.txt \\
            -t ${task.cpus} \\
            -o ${output} \\
            -a ${output_ws}

        get_nextflow_log.bash create_panel_of_normals.log
        """
    }
}
// Variant resources like dbSNP or gNOMAD have no samples by default
