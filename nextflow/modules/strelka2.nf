process STRELKA2 {
    ext version: "2.9.10"

    publishDir "$meta.out", mode: "copy", saveAs: params.saveFn
    publishDir "$meta.log", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(normal), path(tumor), path(indices, arity: "2"), path(manta_indels)
    val(reference)
    val(target_intervals) // For exome data, path to target intervals file
    val(module_number)
    //

    output:
    tuple val(meta), path("*.vcf.gz"), emit: variants
    path(out)
    path("*.log")
    //

    script:
    out = Utils.getName(module_number, meta, "StrelkaOut")
    prefix = Utils.getName(module_number, meta)
    check = file("${meta.out}/${out}")
    target_flag = target_intervals != "" ? " --callRegions=${target_intervals} " : ""
    args = task.ext.args.join(" ")
    if (check.exists()) {
        """
        cp -r ${check} .
        ln -sr ${meta.out}/*_Strelka.vcf.gz .
        ln -sr ${meta.log}/strelka.log .
        """
    } else {
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
    }
    //
}
