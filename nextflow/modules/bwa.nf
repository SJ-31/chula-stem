process BWA {
    ext version: "2.2.1"

    label "big_mem"

    input:


    // Align reads in an alt-aware manner https://gatk.broadinstitute.org/hc/en-us/articles/360037498992--How-to-Map-reads-to-a-reference-with-alternate-contigs-like-GRCH38#1
    tuple val(meta), val(reads)
    val reference

    output:
    //
    tuple val(meta), path(out), emit: mapped
    tuple val(meta), path("*.log"), emit: log

    script:
    //
    out = Utl.getName(meta, null, "bam")
    """
    bwa-mem2 mem \\
        -o aligned.sam \\
        -v 3 \\
        -t ${task.cpus} \\
        ${reference} \\
        ${reads[0]} ${reads[1]}

    samtools view -S -b aligned.sam > aligned.bam
    get_nextflow_log.bash bwa.log

    gatk AddOrReplaceReadGroups \\
        -I aligned.bam \\
        -O ${out} \\
        --RGLB ${meta.RGLB} \\
        --RGPL ${meta.RGPL} \\
        --RGPU ${meta.RGPU} \\
        --RGSM ${meta.RGSM}
    """
}
