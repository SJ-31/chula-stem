process PICARD {
    ext version: params.gatk_version

    publishDir "$meta.out", mode: "copy", saveAs: { x -> x ==~ /.*\.log/ ? null : x }
    publishDir "$meta.log", mode: 'copy', pattern: "*.log"

    input:
    tuple val(meta), path(bam)
    val(omics_type)
    val(reference)
    val(target_intervals)
    val(bait_intervals)
    val(gene_annotations_refFlat)
    val(module_number)

    output:
    tuple (val(meta),
	path("${pair_id}_alignment_metrics.txt"),
	path("${pair_id}_insert_metrics.txt"),
	path("${pair_id}_insert_size_histogram.pdf"),
	path("${pair_id}_depth_out.txt"))

    script:
    out = "${module_number}-${meta.filename}-alignment_metrics_Picard.txt"
    check = file(${meta.out}/"${out}")
    if (check.exists()) {
        """
        ln -sr "${meta.out}/${module_number}"-*_metrics_Picard*" .
        ln -sr "${meta.log}/picard.log" .
        """
    } else {
        """
        gatk CollectAlignmentSummaryMetrics -I $bam \\
            -O $out
        """
        if (omics_type == "exome") {
            """
            gatk CollectHsMetrics -I $bam \\
                --BAIT_INTERVALS $bait_intervals \\
                --TARGET_INTERVALS $target_intervals \\
                -O ${module_number}-${meta.filename}-hs_metrics_Picard.txt
            """
        } else if (omics_type == "wgs") {
            """
            gatk CollectWgsMetrics -I $bam \\
                --REFERENCE_SEQUENCE $reference \\
                -O ${module_number}-${meta.filename}-wgs_metrics_Picard.txt
            """
        } else if (omics_type == "rna-seq") {
            """
            gatk CollectRnaSeqMetrics -I $bam \\
                --REF_FLAT $gene_annotations_refFlat \\
                -O ${module_number}-${meta.filename}-rnaseq_metrics_Picard.txt
            """
        } else {
            throw new Exception("'omics_type' must be one of wgs|rna-seq|exome")
        }
        """
        get_nextflow_log.bash picard.log
        """
    }
}
