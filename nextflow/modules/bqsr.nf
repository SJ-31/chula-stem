process BQSR {
    // GATK version 4.6.1.0
    // Recalibrate and also obtain covariate plot
    input:
    tuple val(meta), path(bam)
    val(reference)
    val(known_variants)
    val(module_number)

    output:
    tuple val(meta), path("*_recal.bam")
    path(report)

    script:
    recal = "${module_number}-${meta.id}_recal.bam"
    report = "${module_number}-${meta.id}_AnalyzeCovariates.pdf"
    def check = file("${meta.out}/$recal")
    def check2 = file("${meta.out}/$report")
    if (check.exists() && check2.exists()) {
        """
        cp $check.name .
        cp $check2.name .
        """
    } else {
        """
        gatk BaseRecalibrator \
            -R $reference \
            -I $bam \
            --known-sites ${known_variants} \
            -O recalibration_1.table

        gatk ApplyBQSR \
            -R $reference \
            -I $bam \
            -bqsr-recal-file recalibration_1.table \
            -O $recal

        gatk BaseRecalibrator \
            -R $reference \
            -I $recal \
            --known-sites ${known_variants} \
            -O recalibration_2.table

        gatk AnalyzeCovariates \
            -before recalibration_1.table \
            -after recalibration_2.table \
            -plots $report
        """
    }

}
