process BQSR {
    ext version: "4.6.1.0"

    publishDir "$meta.out", mode: "copy", saveAs: { x -> x ==~ /.*\.log/ ? null : x }
    publishDir "$meta.log", mode: "copy", pattern: "*.log"

    // Recalibrate and also obtain covariate plot
    input:
    tuple val(meta), path(bam)
    val(reference)
    val(known_sites) // Array containing known sites
    val(module_number)

    output:
    tuple val(meta), path("*_recal.bam"), emit: bam
    path(report)
    path(recal_dir)
    path("*.log")

    script:
    recal = "${module_number}-${meta.id}_recal.bam"
    report = "${module_number}-${meta.id}_AnalyzeCovariates.pdf"
    recal_dir = "${module_number}-recalibration_tables"
    check = file("${meta.out}/$recal")
    check2 = file("${meta.out}/$report")
    def sites_command = known_sites.collect{"--known-sites $it"}.join(' ')
    if (check.exists() && check2.exists()) {
        """
        ln -sr $check .
        ln -sr $check2 .
        ln -sr "${meta.log}/bqsr.log" .
        ln -sr "${meta.out}/${recal_dir}" .
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

        mkdir "${recal_dir}"
        mv recalibration_1.table $recal_dir
        mv recalibration_2.table $recal_dir
        get_nextflow_log.bash bqsr.log
        """
    }

}
