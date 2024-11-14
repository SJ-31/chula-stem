process FASTP {
    ext version: "0.23.4"
    
    publishDir "$meta.out", mode: 'copy'
    publishDir "$meta.log", mode: 'copy', pattern: "*.log"

    input:
    tuple val(meta), val(reads)
    val(module_number) // Helper for tracking module output
    //

    output:
    tuple val(meta), path("*fastp.fastq.gz"), emit: passed
    path("*fail.fastq.gz"), optional: true, emit: failed_reads
    path("*.html")
    path("*.json")
    //

    script:
    n0 = reads[0].baseName.replace(".fastq.gz", "")
    n1 = reads[1].baseName.replace(".fastq.gz", "")
    output1 = "${module_number}-${n1}.fastp.fastq.gz"
    output2 = "${module_number}-${n2}.fastp.fastq.gz"
    prefix = "${module_number}-${meta.id}_fastp"
    check1 = file("${meta.out}/$output1")
    check2 = file("${meta.out}/$output2")
    if (check1.exists() && check2.exists()) {
        """
      ln -sr $check1 .
      ln -sr $check2 .
      ln -sr "${meta.out}/${prefix}.json" .
      ln -sr "${meta.out}/${prefix}.html" .
      ln -sr "${meta.log}/fastp.log" .
        """
    } else {
        """
        fastp -i ${reads[0]} -I ${reads[1]} \\
            -h ${prefix}.html \\
            -j ${prefix}.json \\
            -R ${meta.id}_report \\
            --failed_out ${meta.id}.fail.fastq.gz \\
            -o $output1 -O $output2

        cp .command.out fastp.log
        """
    }

    //
}
