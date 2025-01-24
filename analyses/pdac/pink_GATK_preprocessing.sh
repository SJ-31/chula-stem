GATK="java -Xmx8G -jar /data2/pink/tools/anaconda3/envs/gatk/gatk-4.2.3.0/gatk-package-4.2.3.0-local.jar"
REF="/data2/pink/resource/Homo_sapiens_assembly38.fasta"
## external references
SNP="/data2/pink/resource/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
INDEL="/data2/pink/resource/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
DBSNP="/data2/pink/resource/Homo_sapiens_assembly38.dbsnp138.vcf"
for i in 68
do

        # sample names
        N="N$i"
        NORMAL_SAMPLE="PDAC_${N}_B"

        ## sequencing data files
        WD="/data2/pink/PDAC/PDAC_${N}"
	mkdir ${WD}
	cd ${WD}
        NORMAL1="/mnt/nas/pink/Exome_PDAC/D_P${i}_B/D_P${i}_B_1.trim.fastq.gz"
        NORMAL2="/mnt/nas/pink/Exome_PDAC/D_P${i}_B/D_P${i}_B_2.trim.fastq.gz"

/data2/pink/tools/anaconda3/envs/bwa/bin/bwa mem -t 4 \
        -R "@RG\tID:${N}_EXO_B\tSM:${NORMAL_SAMPLE}\tPL:illumina\tLB:lib_${N}_B\tPU:${N}_B" \
        -M ${REF} ${NORMAL1} ${NORMAL2} > aligned_${NORMAL_SAMPLE}.sam

	${GATK} SortSam \
        -I aligned_${NORMAL_SAMPLE}.sam \
        -O aligned_${NORMAL_SAMPLE}_sorted.bam \
        -SO coordinate

        ${GATK} MarkDuplicates \
          -I aligned_${NORMAL_SAMPLE}_sorted.bam \
          -M aligned_${NORMAL_SAMPLE}_dedup_metrics.txt \
          -O ${NORMAL_SAMPLE}_dedup.bam
        ${GATK} BaseRecalibrator \
           -I ${NORMAL_SAMPLE}_dedup.bam \
           -R ${REF} \
           --known-sites ${SNP} \
           --known-sites ${INDEL} \
           --known-sites ${DBSNP} \
            -O recal_data_${NORMAL_SAMPLE}.table

         ${GATK} ApplyBQSR \
           -R ${REF} \
           -I ${NORMAL_SAMPLE}_dedup.bam \
           --bqsr-recal-file recal_data_${NORMAL_SAMPLE}.table \
           -O recal_data_${NORMAL_SAMPLE}.bam
done

