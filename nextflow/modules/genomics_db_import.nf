process GENOMICS_DB_IMPORT {
    ext version: params.gatk_version

    publishDir "${meta.out}", mode:"copy", saveAs: params.saveFn
    publishDir "${meta.log}", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(normals)
    val(reference)
    val(target_intervals) // In interval list format
    val(module_number)
    //

    output:
    tuple val(meta), path(output), emit: db
    path("*.log")
    //

    script:
    output = Utils.getName(module_number, meta, "Genomics_DB")
    check = file("${meta.out}/${output}")
    normal_flag = normals.collect({" -V ${it} "}).join(" ")
    args = task.ext.args.join(" ")
    if (check.exists()) {
        """
        cp -r ${check} .
        ln -sr ${meta.log}/genomics_db_import.log .
        """
    } else {
        """
        gatk GenomicsDBImport ${args} \\
            -R ${reference} -L ${target_intervals} \\
            --genomicsdb-workspace-path ${output} \\
            ${normal_flag}

        get_nextflow_log.bash genomics_db_import.log
        """
    }
    //
}
