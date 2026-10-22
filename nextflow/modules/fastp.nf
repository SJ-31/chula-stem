process FASTP {
    ext version: "0.23.4"



    errorStrategy "ignore" // BUG: [2025-07-24 Thu] a weird result for one sample causes

    input:
    // this to fail
    tuple val(meta), val(reads)

    output:

    //
    tuple val(meta), path("*fastp.fastq.gz"), emit: passed
    tuple val(meta), path ("*fail.fastq.gz"), optional: true, emit: failed_reads
    path "*.html"
    tuple val(meta), path ("*.json"), emit: json
    tuple val(meta), path("*.log"), emit: logfile

    script:
    //
    n0 = reads[0].baseName.replaceAll(".fastq", "")
    n1 = reads[1].baseName.replaceAll(".fastq", "")
    output1 = "${n0}.fastp.fastq.gz"
    output2 = "${n1}.fastp.fastq.gz"
    prefix = Utl.getName(meta, "fastp")
    args = task.ext.args.join(" ")

    """
    fastp -i "${reads[0]}" -I "${reads[1]}" \\
        -h "${prefix}".html \\
        -j "${prefix}".json \\
        -R "${meta.id}"_report \\
        --failed_out "${meta.id}".fail.fastq.gz \\
        ${args} \\
        -o "${output1}" -O "${output2}"

    get_nextflow_log.bash fastp.log
    """
    //
}
