// Prepare input data for cnv kit
process CNVKIT_PREP {
    ext version: "0.9.11"

    publishDir "$meta.out", mode: "copy", saveAs: { x -> x ==~ /.*\.log/ ? null : x }
    publishDir "$meta.log", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(cohort) // To analyze cohort sequenced on single platform, documentation recommends
    // combining all normal samples into a pooled reference, regardless of whether or not
    // matching tumor-normal pairs were sequenced
    val(reference)
    val(bait_intervals) // bait bed file for exome platforms
    val(exclude) // Genomic regions to exclude
    val(module_number)
    //

    output:
    path(out), emit: reference
    tuple path(target), path(antitarget)
    path(access)
    path("*.log")
    //

    script:
    prefix = "${module_number}-${meta.filename}-Cnvkit"
    out = "${prefix}_reference.cnn"
    target = "${prefix}_target.bed"
    antitarget = "${prefix}_antitarget.bed"
    check = file("${meta.out}/${out}")
    access = "${prefix}_access.bed"
    if (check.exists()) {
        """
        ln -sr $check .
        ln -sr ${meta.out}/${target} .
        ln -sr ${meta.out}/${antitarget} .
        ln -sr ${meta.out}/${access} .
        ln -sr ${meta.log}/cnvkit_prep.log .
        """
    } else {
        """
        cnvkit.py access ${reference} -x ${exclude} -o ${access}
        cnvkit.py target --split ${bait_intervals} -o ${target}
        cnvkit.py antitarget ${target} -g ${access} -o ${antitarget}

        cnvkit.py batch \\
            --normal *.bam \\
            --fasta ${reference} \\
            --output-reference ${out} \\
            --access ${access} \\
            --targets ${target} \\
            --antitargets ${antitarget}

        get_nextflow_log.bash cnvkit_prep.log
        """
    }
    //
}
