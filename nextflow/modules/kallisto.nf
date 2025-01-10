process KALLISTO {
    ext version: "0.51.1"

    publishDir "${meta.out}", mode:"copy", saveAs: params.saveFn
    publishDir "${meta.log}", mode: "copy", pattern: "*.{log,json}"

    input:
    tuple val(meta), val(reads)
    val(index) // Kallisto index of transcriptome fasta file
    val(strandedness) // foward, reverse or null
    val(module_number)
    //

    output:
    tuple val(meta), path(tsv), emit: tsv
    tuple val(meta), path(h5), emit: h5
    path("*.json")
    path("*.log")
    //

    script:
    tsv = Utl.getName(module_number, meta, "kallisto", "tsv")
    h5 = Utl.getName(module_number, meta, "kallisto", "h5")
    check1 = file("${meta.out}/${tsv}")
    check2 = file("${meta.out}/${h5}")
    if (strandedness == "forward") {
        strandedness_flag = " --fr-stranded "
    } else if (strandedness == "reverse") {
        strandedness_flag = " --rf-stranded "
    } else {
        strandedness_flag = ""
    }
    if (check1.exists() && check2.exists()) {
        """
        ln -sr ${check1} .
        ln -sr ${check2} .

        ln -sr ${meta.log}/kallisto_run_info.json .
        ln -sr ${meta.log}/kallisto.log .
        """
    } else {
        """
        kallisto quant --index=${index} \\
            ${strandedness_flag} \\
            --threads=${task.cpus} \\
            --output-dir="." \\
            ${reads[0]} ${reads[1]}

        mv abundance.tsv ${tsv}
        mv abundance.h5 ${h5}

        mv run_info.json kallisto_run_info.json
        get_nextflow_log.bash kallisto.log
        """
    }
    //
}
