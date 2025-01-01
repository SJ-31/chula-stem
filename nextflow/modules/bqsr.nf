process BQSR {
    ext version: params.gatk_version

    label "big_mem"
    publishDir "$meta.out", mode: "copy", saveAs: params.saveFn
    publishDir "$meta.log", mode: "copy", pattern: "*.log"

    // Recalibrate and also obtain covariate plot
    input:
    tuple val(meta), path(bam)
    val(reference)
    val(known_sites) // Array containing known sites
    val(record_recal)
    val(module_number)

    output:
    tuple val(meta), path(recal), emit: bam
    path(report), optional: true
    path(recal_dir)
    path("*.log")

    script:
    recal = Utl.getName(module_number, meta, "recal", "bam")
    report = Utl.getName(module_number, meta, "AnalyzeCovariates", "pdf")
    recal_dir = "${module_number}-recalibration_tables"
    check = file("${meta.out}/$recal")
    check2 = file("${meta.out}/$report")
    def sites_command = known_sites.collect{"--known-sites $it"}.join(' ')
    if (check.exists() && check2.exists()) {
        """
        ln -sr $check .
        ln -sr $check2 .
        ln -sr "${meta.log}/bqsr.log" .
        cp -r "${meta.out}/${recal_dir}" .
        """
    } else {
        """
        gatk BaseRecalibrator \\
            -R $reference \\
            -I $bam \\
            $sites_command \\
            -O recalibration_1.table

        gatk ApplyBQSR \\
            -R $reference \\
            -I $bam \\
            -bqsr-recal-file recalibration_1.table \\
            -O $recal

        mkdir "${recal_dir}"
        mv recalibration_1.table $recal_dir

        if [[ ${record_recal} == "true" ]]; then
            gatk BaseRecalibrator \\
                -R $reference \\
                -I $recal \\
                $sites_command \\
                -O recalibration_2.table

            gatk AnalyzeCovariates \\
                -before recalibration_1.table \\
                -after recalibration_2.table \\
                -plots $report

            mv recalibration_2.table $recal_dir
        else
            echo "Report not requested" > ${report}
        fi

        get_nextflow_log.bash bqsr.log
        """
    }

}
