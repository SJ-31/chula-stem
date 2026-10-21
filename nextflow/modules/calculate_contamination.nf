process CALCULATE_CONTAMINATION {
    ext version: params.gatk_version

    input:
    tuple val(meta), path(pileup)

    output:

    //
    tuple val(meta), path(contamination), emit: c
    tuple val(meta), path(segmentation), emit: s
    tuple val(meta), path("*.log"), emit: log

    script:
    //
    segmentation = Utl.getName(meta, "Segmentation", "table")
    contamination = Utl.getName(meta, "Contamination", "table")
    """
    gatk CalculateContamination \\
        -I ${pileup} \\
        -tumor-segmentation ${segmentation} \\
        -O ${contamination}

    get_nextflow_log.bash calculate_contamination.log
    """
}
