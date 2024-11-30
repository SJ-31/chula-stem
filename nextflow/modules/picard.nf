process PICARD {
    ext version: params.gatk_version

    publishDir "$meta.out", mode: "copy", saveAs: params.saveFn
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
    out = params.getName(module_number, meta, "Picard_alignment_metrics", "txt")
    exome = params.getName(module_number, meta, "Picard_hs_metrics", "txt")
    wgs = params.getName(module_number, meta, "Picard_wgs_metrics", "txt")
    rnaseq = params.getName(module_number, meta, "Picard_rnaseq_metrics", "txt")
    check = file(${meta.out}/"${out}")
    if (check.exists()) {
        """
        ln -sr "${meta.out}/${module_number}"-Picard_*_metrics.txt" .
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
                -O ${exome}
            """
        } else if (omics_type == "wgs") {
            """
            gatk CollectWgsMetrics -I $bam \\
                --REFERENCE_SEQUENCE $reference \\
                -O ${wgs}
            """
        } else if (omics_type == "rna-seq") {
            """
            gatk CollectRnaSeqMetrics -I $bam \\
                --REF_FLAT $gene_annotations_refFlat \\
                -O ${rnaseq}
            """
        } else {
            throw new Exception("'omics_type' must be one of wgs|rna-seq|exome")
        }
        """
        get_nextflow_log.bash picard.log
        """
    }
}
