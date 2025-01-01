process LOFREQ {
    ext version: params.gatk_version

    publishDir "$meta.out", mode: "copy", saveAs: params.saveFn
    publishDir "$meta.log", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(normal), path(tumor), path(indices, arity: "2")
    val(reference)
    val(target_intervals) // For exome data, path to target intervals file
    val(module_number)
    //

    output:
    tuple val(meta), path(out), emit: variants
    path("*.log")
    //

    script:
    out = Utl.getName(module_number, meta, "Mutect2", "vcf.gz")
    raw = Utl.getName(module_number, meta, "Mutect2_raw", "tar.gz")
    check = file("${meta.out}/${out}")
    if (check.exists()) {
        """
        ln -sr $check .
        """
    } else {
        """
        lofreq somatic -n ${normal} -t ${tumor} -f ${reference} \\
            --threads 8 -o out_ \\
            -d ${known_variants}

        get_nextflow_log.bash lofreq.log
        """
    }
    //
}
