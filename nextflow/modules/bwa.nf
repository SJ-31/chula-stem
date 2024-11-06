process BWA {
    // bwa-mem2 version 2.2.1
    publishDir "$meta.out", mode: "copy"
    publishDir "$meta.log", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), val(reads)
    val(reference)
    val(module_number)
    //

    output:
    tuple val(meta), path(out), emit: mapped
    path("*.log")
    //

    script:
    out = "${module_number}-${meta.id}.bam"
    def check = file("${meta.out}/$out")
    if (check.exists()) {
        """
        cp $check.name .
        cp "${meta.log}/bwa.log" .
        """
    } else {
        """
        bwa-mem2 index $reference
        bwa-mem2 mem \
            -o aligned.sam \
            -v 3 \
            $reference \
            ${reads[0]} ${reads[1]}

        samtools view -S -b aligned.sam > $out
        cp .command.out bwa.log
        """
    }
    //
}
