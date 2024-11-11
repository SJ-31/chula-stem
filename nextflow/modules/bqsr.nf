process BQSR {
    ext version: "4.6.1.0"
    conda { task.ext.conda }
    // Recalibrate and also obtain covariate plot
    input:
    tuple val(meta), path(bam)
    val(reference)
    val(known_sites) // Array containing known sites
    val(module_number)

    output:
    tuple val(meta), path("*_recal.bam"), emit: bam
    path(report)

    script:
    recal = "${module_number}-${meta.id}_recal.bam"
    report = "${module_number}-${meta.id}_AnalyzeCovariates.pdf"
    check = file("${meta.out}/$recal")
    check2 = file("${meta.out}/$report")
    def sites_command = known_sites.collect{"--known-sites $it"}.join(' ')
    if (check.exists() && check2.exists()) {
        """
        cp $check .
        cp $check2 .
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

        gatk BaseRecalibrator \\
            -R $reference \\
            -I $recal \\
            $sites_command \\
            -O recalibration_2.table

        gatk AnalyzeCovariates \\
            -before recalibration_1.table \\
            -after recalibration_2.table \\
            -plots $report
        """
    }

}
