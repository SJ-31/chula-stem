process STAR_FUSION {
    ext version: "1.13.0"

    input:
    tuple val(meta), path(chimeric_junction)
    val genome_lib // Library containing fusion annotation info

    output:

    //
    tuple val(meta), path(fusion)
    path abridged
    tuple val(meta), path("*.log"), emit: logfile

    script:
    //
    fusion = Utl.getName(meta, "Star_Fusion", "tsv")
    abridged = Utl.getName(meta, "Star_Fusion_abridged", "tsv")
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
