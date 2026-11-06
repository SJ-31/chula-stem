process STRELKA2 {
    ext version: "2.9.10"

    input:
    tuple val(meta), path(normal), path(tumor), path(indices, arity: "2"), path(manta_indels)
    val reference
    val target_intervals // For exome data, path to target intervals file

    output:

    //
    tuple val(meta), path("*.vcf.gz"), emit: variants
    tuple val(meta), path(out), emit: dir
    tuple val(meta), path("*.log"), emit: logfile

    script:
    //
    out = Utl.getName(meta, "StrelkaOut")
    prefix = Utl.getName(meta)
    target_flag = target_intervals != "" ? " --callRegions=${target_intervals} " : ""
    args = task.ext.args.join(" ")
    """
    tabix -f ${manta_indels}

    configureStrelkaSomaticWorkflow.py \\
        --normalBam ${normal} \\
        --tumorBam ${tumor} \\
        --referenceFasta ${reference} \\
        --indelCandidates ${manta_indels} \\
        ${target_flag} \\
        ${args} \\
        --runDir ${out}

    ${out}/runWorkflow.py -m local

    rename_strelka.bash -n "${meta.RGSM_normal}" -t "${meta.RGSM_tumor}" \\
        -o "${out}" -e "${params.source_name}" -d "${params.source_description}" -p "${prefix}"

    get_nextflow_log.bash strelka.log
    """
    //
}
