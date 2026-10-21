process CALCULATE_CONTAMINATION {
    ext version: params.gatk_version

    
    

    input:
    tuple val(meta), path(pileup)
    
    //

    output:
    tuple val(meta), path(contamination), emit: c
    tuple val(meta), path(segmentation), emit: s
    tuple val(meta), path("*.log"), emit: log
    //

    script:
    segmentation = Utl.getName(meta, "Segmentation", "table")
    contamination = Utl.getName(meta, "Contamination", "table")
    check1 = file("${meta.out}/${segmentation}")
    check2 = file("${meta.out}/${contamination}")
    args = task.ext.args.join(" ")
    if (check1.exists() && check2.exists()) {
        """
        ln -sr ${check1} .
        ln -sr ${check2} .
        ln -sr ${meta.log}/calculate_contamination.log .
        """
    } else {
        """
        gatk CalculateContamination \\
            -I ${pileup} \\
            -tumor-segmentation ${segmentation} \\
            -O ${contamination}

        get_nextflow_log.bash calculate_contamination.log
        """
    }
    //
}
