process BQSR {
    ext version: params.gatk_version

    label "big_mem"
    
    

    // Recalibrate and also obtain covariate plot
    input:
    tuple val(meta), path(bam)
    val(reference)
    val(known_sites) // Array containing known sites
    val(record_recal)
    

    output:
    tuple val(meta), path(recal), emit: bam
    path(report), emit: report, optional: true
    tuple val(meta), path(recal_dir), emit: recal_tables
    tuple val(meta), path("*.log"), emit: log

    script:
    recal = Utl.getName(meta, "recal", "bam")
    report = Utl.getName(meta, "AnalyzeCovariates", "pdf")
    recal_dir = "recalibration_tables"
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
