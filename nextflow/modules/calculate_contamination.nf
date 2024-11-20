process CALCULATE_CONTAMINATION {
    ext version: params.gatk_version

    publishDir "$meta.out", mode: "copy", saveAs: { x -> x ==~ /.*\.log/ ? null : x }
    publishDir "$meta.log", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(pileup)
    val(module_number)
    //

    output:
    tuple val(meta), path(contamination), emit: c
    tuple val(meta), path(segmentation), emit: s
    path("*.log")
    //

    script:
    segmentation = "${module_number}-${meta.filename}-Segmentation.table"
    contamination = "${module_number}-${meta.filename}-Contamination.table"
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
