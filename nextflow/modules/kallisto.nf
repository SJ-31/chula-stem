process KALLISTO {
    ext version: "0.51.1"

    input:
    tuple val(meta), val(reads)
    val index // Kallisto index of transcriptome fasta file
    val strandedness // foward, reverse or null

    output:

    //
    tuple val(meta), path(tsv), emit: tsv
    tuple val(meta), path(h5), emit: h5
    path "*.json"
    tuple val(meta), path("*.log"), emit: log

    script:
    //
    tsv = Utl.getName(meta, "kallisto", "tsv")
    h5 = Utl.getName(meta, "kallisto", "h5")
    if (strandedness == "forward") {
        strandedness_flag = " --fr-stranded "
    }
    else if (strandedness == "reverse") {
        strandedness_flag = " --rf-stranded "
    }
    else {
        strandedness_flag = ""
    }

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
