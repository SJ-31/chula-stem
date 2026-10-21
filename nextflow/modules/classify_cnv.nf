process CLASSIFY_CNV {
    ext version: "1.1.1"

    input:
    tuple val(meta), path(bed)

    output:

    //
    // Input format for the BED file is to have the columns below
    // 1. chromosome
    // 2. CNV start
    // 3. CNV end
    // 4. CNV type (DEL or DUP)
    tuple val(meta), path(tsv), emit: tsv
    path intermediate
    tuple val(meta), path("*.log"), emit: log

    script:
    //
    prefix = Utl.getName(meta, "ClassifyCNV")
    tsv = "${prefix}.tsv"
    intermediate = "${prefix}_intermediate"
    """
    ClassifyCNV.py --infile ${bed} \\
        --GenomeBuild ${task.ext.genome} \\
        --outdir ${intermediate}

    mv "${intermediate}/Scoresheet.txt" "${prefix}.tsv"
    get_nextflow_log.bash classify_cnv.log
    """
}
