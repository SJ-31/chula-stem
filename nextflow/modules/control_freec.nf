process CONTROL_FREEC {
    ext version: "11.6b"

    // ABANDONED
    label "mid_mem"
    publishDir "${meta.out}", mode:"copy", saveAs: params.saveFn
    publishDir "${meta.log}", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(normal), path(tumor), path(indices, arity: "2")
    val(reference)
    val(known_snps) // Specifically asks for file with known snps
    val(targets) // Target interval bed file for exome data
    val(module_number)
    //

    output:
    tuple val(meta), path(cnvs), emit: cnvs
    path(prefix)
    path("*.log")
    //

    script:
    // See https://boevalab.inf.ethz.ch/FREEC/tutorial.html#OUTPUT to interpret output
    tname = tumor.name
    nname = normal.name

    prefix = "${module_number}-${meta.filename}-Freec"
    info = "${prefix}_info.txt"
    cnvs = "${prefix}_CNVs"
    ratios = "${prefix}_ratio.txt"
    tumor_cpn = "${prefix}_tumor.cpn"
    normal_cpn = "${prefix}_normal.cpn"

    check = file("${meta.out}/${prefix}")
    interval_flag = targets != "" ? "--intervals ${targets}" : ""
    args = task.ext.args.join(" ")
    if (check.exists()) {
        """
        cp -r ${check} .
        ln -sr ${meta.log}/freec.log .
        """
    } else {
        """
        split_contigs.bash ${reference} chrs

        make_freec_config --chr_len_file "${reference}.fai" \\
            --read_type paired_end \\
            --normal ${normal} \\
            --tumor ${normal} \\
            --file_format pileup \\
            --chr_files chrs \\
            --sex ${meta.sex} \\
            --snps ${known_snps} \\
            ${interval_flag} \\
            --output config \
            --noisy_data

        freec -conf config

        mkdir ${prefix}
        mv "${tname}_info.txt" "${prefix}/${info}"
        mv "${tname}_CNVs" "${prefix}/${cnvs}"
        mv "${tname}_ratio.txt" "${prefix}/${ratios}"
        mv "${tname}_sample.cpn" ${prefix}/${tumor_cpn}
        mv "${nname}_control.cpn" "${prefix}/${normal_cpn}"

        get_nextflow_log.bash freec.log
        """
    }
    //
}
