process CLASSIFY_CNV {
    ext version: "1.1.1"

    publishDir "${meta.out}", mode:"copy", saveAs: params.saveFn
    publishDir "${meta.log}", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(bed)
    val(module_number)
    //
    // Input format for the BED file is to have the columns below
    // 1. chromosome
    // 2. CNV start
    // 3. CNV end
    // 4. CNV type (DEL or DUP)

    output:
    tuple val(meta), path(tsv), emit: tsv
    path(intermediate)
    path("*.log")
    //

    script:
    prefix = Utl.getName(module_number, meta, "ClassifyCNV")
    tsv = "${prefix}.tsv"
    intermediate = "${prefix}_intermediate"
    check = file("${meta.out}/${tsv}")
    args = task.ext.args.join(" ")
    if (check.exists()) {
        """
        ln -sr ${check} .
        cp -r "${meta.out}/${intermediate}" .
        ln -sr ${meta.log}/classify_cnv.log .
        """
    } else {
        """
        ClassifyCNV.py --infile ${bed} \\
            --GenomeBuild ${task.ext.genome} \\
            --outdir ${intermediate}

        mv "${intermediate}/Scoresheet.txt" "${prefix}.tsv"
        get_nextflow_log.bash classify_cnv.log
        """
    }
    //
}
