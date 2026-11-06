nextflow.enable.types = true

process MERGE_FASTQ {
    input:
    record(id: String,
           fastq_dir: String,
           status: String
    )

    output:
    record(id: id,
           status: status,
           fastq_1: file(fq1),
           fastq_2: file(fq1)
    )
    script:
    fq1 = "${id}_R1.fastq.gz"
    fq2 = "${id}_R2.fastq.gz"
    """
    seqkit ${fastq_dir}/*_R1.*fastq.gz > t1.fq.gz
    mv t1.fq.gz ${fq1}
    
    seqkit ${fastq_dir}/*_R2.*fastq.gz > t2.fq.gz
    mv t2.fq.gz ${fq2}
    """
}
