nextflow.enable.types = true

process MERGE_FASTQ {
    input:
    record(
        id: String,
        fastq_dir: String,
        status: Integer
    )

    output:
    record(
        id: id,
        status: status,
        fastq_1: file(fq1),
        fastq_2: file(fq2),
    )

    script:
    fq1 = "${id}_R1.fastq.gz"
    fq2 = "${id}_R2.fastq.gz"
    """
    seqkit seq ${fastq_dir}/*_R1.*fastq.gz -o t1.fq.gz
    mv t1.fq.gz ${fq1}
    
    seqkit seq ${fastq_dir}/*_R2.*fastq.gz -o t2.fq.gz
    mv t2.fq.gz ${fq2}
    """
}
