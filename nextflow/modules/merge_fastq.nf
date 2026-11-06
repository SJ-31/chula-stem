nextflow.enable.types = true

process MERGE_FASTQ {
    input:
    record(id: String,
           fastq_dir: Path,
           status: String
    )

    output:
    record(id: id,
           status: status,
           fastq_1: fq1,
           fastq_2: fq2
    )
    script:
    fq1 = "${id}_R1.fastq.gz"
    fq2 = "${id}_R2.fastq.gz"
    """
    seqkit *_R1.*fastq.gz > t1.fq.gz
    mv t1.fq.gz ${fq1}
    
    seqkit *_R2.*fastq.gz > t2.fq.gz
    mv t2.fq.gz ${fq2}
    """
}
