process FACETS {
    ext version: "2"

    publishDir "$meta.out", mode: "copy", saveAs: params.saveFn
    publishDir "$meta.log", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(counts_file) // Output from facets_pileup
    val(cnvkit_autobin) // estimate of bin size produced by `cnvkit.py autobin` stdout
    val(module_number)
    //

    output:
    tuple val(meta), path(prefix), emit: cnv
    tuple val(with_caller), path("${prefix}/${prefix}_hisens.rds"), emit: rds
    tuple val(meta), path("${prefix}/purity.txt"), path("${prefix}/ploidy.txt"), emit: purity_ploidy
    path(tsv)
    path("*.log")
    //

    script:
    prefix = Utl.getName(module_number, meta, "Facets")
    check = file("${meta.out}/${prefix}")
    tsv = "${prefix}/${prefix}_hisens.tsv"
    check2 = file("${meta.out}/${tsv}")
    with_caller = meta + ["caller": "facets"]
    if (cnvkit_autobin) {
        bin_size = file(cnvkit_autobin).readLines()[1].split()[2]
        size_flag = Utl.overrideArgs(["--snp-window-size ${bin_size}"], task.ext.args)
    } else {
        size_flag = ""
    }
    size_flag = ""
    args = task.ext.args.join(" ")

    if (check.exists() && check2.exists()) {
        """
        cp -r $check .
        ln -sr ${check2} .
        ln -sr ${meta.log}/facets.log .
        """
    } else {
        """
        run-facets-wrapper.R \\
            --counts-file ${counts_file} \\
            --sample-id ${prefix} \\
            --genome ${task.ext.genome} \\
            ${size_flag} \\
            ${args}

        cut -f 2 ${prefix}/${prefix}.txt | tail -n 1 > ${prefix}/purity.txt
        cut -f 3 ${prefix}/${prefix}.txt | tail -n 1 > ${prefix}/ploidy.txt

        get_facets.bash "${prefix}/${prefix}_hisens.rds" "${tsv}"

        get_nextflow_log.bash facets.log
        """
    }
    //
}
// <2025-02-27 Thu> There's a bug in the "calculate_ntai" function of copy_number_scores.R
// This is being
// It's happening in this snippet:
//  # Check relative position to centromere
// if(chrom_segs$AI[1] == 2 & nrow(chrom_segs) != 1 & chrom_segs$end[1] < (sample_chrom_info$centromere[chr])){
//     segs$AI[which(segs$chrom == chr)][1] = 1 # if the first segment of chromosome is AI and does not extend to centromere --> telomeric AI
// }
