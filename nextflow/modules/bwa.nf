process BWA {
    ext version: "2.2.1"
   
    label "big_mem"
    publishDir "$meta.out", mode: "copy", saveAs: { x -> x ==~ /.*\.log/ ? null : x }
    publishDir "$meta.log", mode: "copy", pattern: "*.log"
    // Align reads in an alt-aware manner https://gatk.broadinstitute.org/hc/en-us/articles/360037498992--How-to-Map-reads-to-a-reference-with-alternate-contigs-like-GRCH38#1

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
    check = file("${meta.out}/$out")
    if (check.exists()) {
        """
        ln -sr $check .
        ln -sr "${meta.log}/bwa.log" .
        """
    } else {
        """
        bwa-mem2 mem \\
            -o aligned.sam \\
            -v 3 \\
            $reference \\
            ${reads[0]} ${reads[1]}

        samtools view -S -b aligned.sam > aligned.bam
        get_nextflow_log.bash bwa.log

        gatk AddOrReplaceReadGroups \\
            -I aligned.bam \\
            -O $out \\
            --RGLB $meta.RGLB \\
            --RGPL $meta.RGPL \\
            --RGPU $meta.RGPU \\
            --RGSM $meta.RGSM
        """
    }
    //
}
