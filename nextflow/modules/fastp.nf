process FASTP {
    publishDir "$meta.out", mode: 'copy'
    publishDir "$meta.log", mode: 'copy', pattern: "*.log"

    input:
    tuple val(meta), val(reads)
    val(module_number) // Helper for tracking module output
    //

    output:
    tuple val(meta), path("*fastp.fastq.gz"), emit: out
    path("*fail.fastq.gz"), optional: true, emit: failed_reads
    path("*.html")
    path("*.json")
    //

    script:
    output1 = "${reads[0].baseName}.fastp.fastq.gz"
    output2 = "${reads[1].baseName}.fastp.fastq.gz"
    """
    fastp -i ${reads[0]} -I ${reads[1]} \
        -z 4 \
        -h ${module_number}-${meta.sample}_fastp.html \
        -j ${module_number}-${meta.sample}_fastp.json \
        -R ${meta.sample}_report \
        --failed_out ${meta.sample}.fail.fastq.gz \
        -o $output1 -O $output2

    cp .command.out fastp.log
    """
    //
}
