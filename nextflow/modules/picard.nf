process PICARD {
    ext version: params.gatk_version

    publishDir "$meta.out", mode: "copy", saveAs: params.saveFn
    publishDir "$meta.log", mode: 'copy', pattern: "*.log"

    input:
    tuple val(meta), path(bam), path(index)
    val(omics_type) // One of `hs` (exome), `wgs` or `rnaseq`
    val(reference)
    val(target_intervals) // TODO: Requires interval_list format, generated
    val(bait_intervals)
    val(gene_annotations_refFlat)
    val(module_number)

    output:
    tuple path(out), path(out2), emit: metrics
    path("*.log")

    script:
    out = Utils.getName(module_number, meta, "Picard_alignment_metrics", "txt")
    out2 = Utils.getName(module_number, meta, "Picard_${omics_type}_metrics", "txt")
    check = file("${meta.out}/${out}")
    check2 = file("${meta.out}/${out2}")
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
