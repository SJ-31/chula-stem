gatk CollectAlignmentSummaryMetrics -I "!{bam}" \
    -O "!{out}"
if [[ "!{omics_type}" == "hs" ]]; then
    gatk CollectHsMetrics -I "!{bam}" \
        --BAIT_INTERVALS "!{bait_intervals}" \
        --TARGET_INTERVALS "!{target_intervals}" \
        -O "!{exome}"
elif [[ "!{omics_type}" == "wgs" ]]; then
    gatk CollectWgsMetrics -I "!{bam}" \
        --REFERENCE_SEQUENCE "!{reference}" \
        -O "!{wgs}"
elif [[ "!{omics_type}" == "rnaseq" ]]; then
    gatk CollectRnaSeqMetrics -I !{}bam \
        --REF_FLAT "!{gene_annotations_refFlat}" \
        -O "!{rnaseq}"
else
    echo "'omics_type' must be one of wgs|rna-seq|exome"
    exit 1
fi

get_nextflow_log.bash picard.log
