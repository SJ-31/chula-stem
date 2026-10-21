process CREATE_PANEL_OF_NORMALS {
    ext version: "1"

    input:
    tuple val(meta), val(vcfs)
    val minimum
    val reference
    val other_vcfs

    output:

    //
    tuple val(meta), path(output), emit: pon // Empty pon with no sample information
    tuple val(meta), path(output_ws), emit: pon_ws
    path "*.tbi"
    tuple val(meta), path("*.log"), emit: log

    script:
    //
    output = Utl.getName(meta, "PON", "vcf.gz")
    output_ws = Utl.getName(meta, "PON_ws", "vcf.gz")
    to_sample_spec = vcfs.toList().join("\n")
    to_sample_spec = "${to_sample_spec}\n"
    other_flag = other_vcfs ? " -v other.txt " : ""
    to_other_spec = other_vcfs.join("\n")
    """
    echo -e -n "${to_sample_spec}" > samples.txt
    echo -e -n "${to_other_spec}" > other.txt

    create_panel_of_normals.bash -l samples.txt \\
        -m ${minimum} \\
        ${other_flag} \\
        -t ${task.cpus} \\
        -o tmp.vcf.gz \\
        -a tmp_ws.vcf.gz

    standardize_vcf_clean.bash -i tmp.vcf.gz \\
        -r ${reference} \\
        -o ${output} 

    gatk IndexFeatureFile -I ${output}

    standardize_vcf_clean.bash -i tmp_ws.vcf.gz \\
        -r ${reference} \\
        -o ${output_ws}
    
    gatk IndexFeatureFile -I ${output_ws}
    
    get_nextflow_log.bash create_panel_of_normals.log
    """
}
// Variant resources like dbSNP or gNOMAD have no samples by default
