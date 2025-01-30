process PICARD {
    ext version: params.gatk_version

    publishDir "$meta.out", mode: "copy", saveAs: params.saveFn
    publishDir "$meta.log", mode: 'copy', pattern: "*.log"

    input:
    tuple val(meta), path(bam), path(index)
    val(omics_type) // One of `hs` (exome), `wgs` or `rnaseq`
    val(reference)
    val(target_intervals)
    val(bait_intervals)
    val(strandedness)
    val(gene_annotations_refFlat)
    val(module_number)
    // Both baits and targets must be in interval list format

    output:
    tuple path(out), path(out2), emit: metrics
    path("*.log")

    script:
    out = Utl.getName(module_number, meta, "Picard_alignment_metrics", "txt")
    out2 = Utl.getName(module_number, meta, "Picard_${omics_type}_metrics", "txt")
    check = file("${meta.out}/${out}")
    check2 = file("${meta.out}/${out2}")
    if (strandedness == "reverse" ) {
       strandedness_flag = "SECOND_READ_TRANSCRIPTION_STRAND "
    } else if (strandedness == "forward" || strandedness == "unstranded") {
       strandedness_flag = "FIRST_READ_TRANSCRIPTION_STRAND"
    } else {
       strandedness_flag =  ""
    }
    if (check.exists() && check2.exists()) {
        """
        ln -sr ${check} .
        ln -sr ${check2} .
        ln -sr "${meta.log}/picard.log" .
        """
    } else {
        template "picard.bash"
    }
}
