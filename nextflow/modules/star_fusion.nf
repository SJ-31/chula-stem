process STAR_FUSION {
    ext version: "1.13.0"

    publishDir "${meta.out}", mode:"copy", saveAs: params.saveFn
    publishDir "${meta.log}", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(chimeric_junction)
    val(genome_lib) // Library containing fusion annotation info
    val(module_number)
    //

    output:
    tuple val(meta), path(fusion)
    path(abridged)
    path("*.log")
    //

    script:
    fusion = Utl.getName(module_number, meta, "Star_Fusion", "tsv")
    abridged = Utl.getName(module_number, meta, "Star_Fusion_abridged", "tsv")
    check1 = file("${meta.out}/${fusion}")
    check2 = file("${meta.out}/${abridged}")
    if (check1.exists() && check2.exists()) {
        """
        ln -sr ${check1} .
        ln -sr ${check2} .
        ln -sr ${meta.log}/star_fusion.log .
        """
    } else {
        """
        STAR-Fusion --chimeric_junction ${chimeric_junction} \\
            --genome_lib_dir ${genome_lib} \\
            --output_dir . \\
            --show_full_usage_info \\
            --CPU ${task.cpus}

        mv star-fusion.fusion_predictions.tsv ${fusion}
        mv star-fusion.fusion_predictions.abridged.tsv ${abridged}

        get_nextflow_log.bash star_fusion.log
        """
    }
    //
}
