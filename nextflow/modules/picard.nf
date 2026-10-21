process PICARD {
    ext version: params.gatk_version

    input:
    tuple val(meta), path(bam), path(index)
    val omics_type // One of `hs` (exome), `wgs` or `rnaseq`
    val reference
    val target_intervals
    val bait_intervals
    val strandedness
    val gene_annotations_refFlat

    output:

    // Both baits and targets must be in interval list format
    tuple val(meta), path(out2), emit: specific_metrics
    tuple val(meta), path(out), emit: algn_metrics
    tuple path(out), path(out2), emit: metrics
    tuple val(meta), path("*.log"), emit: log

    script:
    out = Utl.getName(meta, "Picard_alignment_metrics", "txt")
    out2 = Utl.getName(meta, "Picard_${omics_type}_metrics", "txt")
    if (strandedness == "reverse") {
        strandedness_flag = "SECOND_READ_TRANSCRIPTION_STRAND"
    }
    else if (strandedness == "forward" || strandedness == "unstranded") {
        strandedness_flag = "FIRST_READ_TRANSCRIPTION_STRAND"
    }
    else {
        strandedness_flag = ""
    }
    if (!strandedness_flag && omics_type == "rnaseq") {
        throw new Exception("Strandedness must be supplied for RNAseq data")
    }
    """
    gatk CollectAlignmentSummaryMetrics -I "${bam}" \\
        -O "${out}"
    if [[ "${omics_type}" == "hs" ]]; then
        gatk CollectHsMetrics -I "${bam}" \\
            --BAIT_INTERVALS "${bait_intervals}" \\
            --TARGET_INTERVALS "${target_intervals}" \\
            -O "${out2}"
    elif [[ "${omics_type}" == "wgs" ]]; then
        gatk CollectWgsMetrics -I "${bam}" \\
            --REFERENCE_SEQUENCE "${reference}" \\
            -O "${out2}"
    elif [[ "${omics_type}" == "rnaseq" ]]; then
        gatk CollectRnaSeqMetrics -I "${bam}" \\
            --STRAND_SPECIFICITY ${strandedness_flag} \\
            --REF_FLAT "${gene_annotations_refFlat}" \\
            -O "${out2}"
    else
        echo "'omics_type' must be one of wgs|rna-seq|exome"
        exit 1
    fi

    get_nextflow_log.bash picard.log
    """
}
