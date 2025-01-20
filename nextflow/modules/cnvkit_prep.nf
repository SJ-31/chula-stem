// Prepare input data for cnv kit
process CNVKIT_PREP {
    ext version: "0.9.11"

    publishDir "$meta.out", mode: "copy", saveAs: params.saveFn
    publishDir "$meta.log", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(cohort, arity: "1..*") // To analyze cohort sequenced on single platform, documentation recommends
    // combining all normal samples into a pooled reference, regardless of whether or not
    // matching tumor-normal pairs were sequenced
    val(reference)
    val(bait_intervals) // bait bed file for exome platforms. The cli flag for this is confusingly called "targets", but see [1]
    val(exclude) // Genomic regions to exclude
    val(with_normals) // If true, reference is constructed from provided normal files
    // in `cohort`. Otherwise, the a flat reference is made
    val(omics_type)  // wes|amplicon|wgs
    val(module_number)
    //

    output:
    path(out), emit: reference
    tuple path(target), path(antitarget)
    path(access)
    path(bins), emit: autobin, optional: true
    path("*.log")
    //

    script:
    prefix = Utl.getName(module_number, meta, "CnvkitPrep")
    out = "${prefix}_reference.cnn"
    target = "${prefix}_target.bed"
    antitarget = "${prefix}_antitarget.bed"
    bins = "${prefix}_bins.tsv"
    check = file("${meta.out}/${out}")
    omics_type = omics_type == "wes" ? "hybrid" : omics_type
    access = "${prefix}_access.bed"
    if (check.exists()) {
        """
        ln -sr $check .
        ln -sr ${meta.out}/${target} .
        ln -sr ${meta.out}/${antitarget} .
        ln -sr ${meta.out}/${access} .
        ln -sr ${meta.log}/cnvkit_prep.log .
        if [[ -e ${meta.out}/${bins} ]]; then
            ln -sr ${meta.out}/${bins} .
        fi
        """
    } else {
        """
        cnvkit.py access ${reference} -x ${exclude} -o ${access}

        if [[ "${with_normals}" == "true" ]]; then
            cnvkit.py autobin *.bam \\
                --fasta ${reference} \\
                --method ${omics_type} \\
                --access ${access} \\
                --targets ${bait_intervals} \\
                --target-output-bed ${target} \\
                --antitarget-output-bed ${antitarget} > ${bins}
        else
            cnvkit.py target --split ${bait_intervals} -o ${target}
            cnvkit.py antitarget ${target} -g ${access} -o ${antitarget}
        fi
        cnvkit.py reference --output ${out} \\
            --fasta ${reference} \\
            --targets ${target} \\
            --antitargets ${antitarget}

        get_nextflow_log.bash cnvkit_prep.log
        """
    }
    //
}

// [1] https://github.com/etal/cnvkit/blob/master/doc/quickstart.rst#build-a-reference-from-normal-samples-and-infer-tumor-copy-ratios
