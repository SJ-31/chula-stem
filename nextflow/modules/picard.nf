process PICARD {
    ext version: "4.6.1.0"
  

    publishDir "$meta.out", mode: 'copy'
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
    out = "${module_number}-${meta.id}_alignment_metrics_Picard.txt"
    check = file(${meta.out}/"${out}")
    if (check.exists()) {
        """
        cp "${meta.out}/${module_number}"-*_metrics_Picard*" .
        cp "${meta.log}/picard.log" .
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
                -O ${module_number}-${meta.id}_hs_metrics_Picard.txt
            """
        } else if (omics_type == "wgs") {
            """
            gatk CollectWgsMetrics -I $bam \\
                --REFERENCE_SEQUENCE $reference \\
                -O ${module_number}-${meta.id}_wgs_metrics_Picard.txt
            """
        } else if (omics_type == "rna-seq") {
            """
            gatk CollectRnaSeqMetrics -I $bam \\
                --REF_FLAT $gene_annotations_refFlat \\
                -O ${module_number}-${meta.id}_rnaseq_metrics_Picard.txt
            """
        } else {
            throw new Exception("'omics_type' must be one of wgs|rna-seq|exome")
        }
        """
        cp .command.out picard.log
        """
    }
}
