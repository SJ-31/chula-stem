for i in 63_2 #79_1 64_2 70_2 #73_1 62_2 78_1 _63_2 #61_2
do
REF="/data2/pink/resource/Homo_sapiens_assembly38.fasta"
GATK="/data2/pink/tools/anaconda3/envs/gatk/gatk-4.2.3.0/gatk-package-4.2.3.0-local.jar"
## external references
VAR="/data2/pink/resource/af-only-gnomad.hg38.vcf.gz"
PON="/data2/pink/Exome_CEN2/PON/PON2/pon_hg38_edit5.vcf.gz"

## interval for GetPileupSummaries -> CalculateContaminant
INTERVAL_LIST="/data2/pink/Exome_CEN2/S07604514_Covered.interval_list"
TMP="/mnt/nas/pink/tmp"
# sample names
N="N${i}"
NORMAL_SAMPLE="${N}_B"
TUMOR_SAMPLE="PDAC_${N}_C"
DIR="/mnt/nas/pink/result_Exome_HCC/${N}"
DIRR="/data2/pink/PDAC/${N}"
#OD="/data2/pink/Exome_CEN2/EXOME/${N}"
cd ${DIR}
export _JAVA_OPTIONS=-Djava.io.tmpdir=/data/tmp
#java -Xmx30G -jar ${GATK} Mutect2 \
   -R ${REF} \
   -I ${DIRR}/recal_data_${NORMAL_SAMPLE}.bam \
   -I ${DIR}/recal_data_${TUMOR_SAMPLE}.bam \
   -normal ${NORMAL_SAMPLE} \
   --germline-resource ${VAR} \
   --tmp-dir ${TMP} \
   --panel-of-normals ${PON} \
   -L ${INTERVAL_LIST} \
   -O ${DIR}/somatic_variants_${TUMOR_SAMPLE}.vcf.gz

#java -Xmx8G -jar ${GATK} GetPileupSummaries \
    -I ${DIRR}/recal_data_${NORMAL_SAMPLE}.bam \
    -V ${VAR} \
    --tmp-dir ${TMP} \
    -L ${INTERVAL_LIST} \
    -O ${DIR}/${NORMAL_SAMPLE}_pileups.table

#java -Xmx8G -jar ${GATK} GetPileupSummaries \
    -I ${DIR}/recal_data_${TUMOR_SAMPLE}.bam \
    -V ${VAR} \
    --tmp-dir ${TMP} \
    -L ${INTERVAL_LIST} \
    -O ${DIR}/${TUMOR_SAMPLE}_pileups.table

#java -Xmx8G -jar ${GATK} CalculateContamination \
   -I ${DIR}/${TUMOR_SAMPLE}_pileups.table \
   --matched-normal ${DIR}/${NORMAL_SAMPLE}_pileups.table \
   -segments ${DIR}/PDAC_${N}_segments.tsv \
   --tmp-dir ${TMP} \
   -O ${DIR}/PDAC_${N}_contamination.table

#java -Xmx8G -jar ${GATK} FilterMutectCalls \
  -R ${REF} \
  -V ${DIR}/somatic_variants_${TUMOR_SAMPLE}.vcf.gz \
  --contamination-table ${DIR}/PDAC_${N}_contamination.table \
  --tmp-dir ${TMP} \
  --tumor-segmentation ${DIR}/PDAC_${N}_segments.tsv \
  -O ${DIR}/${TUMOR_SAMPLE}_filtered_variants.vcf.gz

#java -Xmx8G -jar ${GATK} SelectVariants \
  -R ${REF} \
  -V ${DIR}/${TUMOR_SAMPLE}_filtered_variants.vcf.gz \
  -O ${DIR}/${TUMOR_SAMPLE}_EXO_remove_filtered.vcf.gz \
  --tmp-dir ${TMP} \
  --create-output-variant-index true \
  --exclude-filtered true \
--exclude-non-variants true

java -Xmx8G -jar ${GATK} Funcotator \
  -V ${DIR}/${TUMOR_SAMPLE}_EXO_remove_filtered.vcf.gz \
  -R ${REF} \
  -O ${DIR}/${TUMOR_SAMPLE}_filtered_variants_funco_selected.maf \
   --output-file-format MAF \
   --ref-version hg38 \
   --tmp-dir ${TMP} \
   --data-sources-path /data2/pink/resource/funcotator_dataSources.v1.6.20190124s \
   --transcript-selection-mode ALL \
   --annotation-default tumor_barcode:${TUMOR_SAMPLE}

#java -Xmx8G -jar ${GATK} LiftoverVcf \
  #-I ${DIR}/${TUMOR_SAMPLE}_EXO_remove_filtered.vcf.gz \
  #-O ${DIR}/PDAC_${N}_EXO_remove_filtered.lift0ver.vcf.gz \
  #-C /data2/pink/Exome_CEN2/EXOME/command/hg38ToHg19.over.chain.gz \
  #--REJECT rejected_variants_PDAC_${N}_rm.vcf \
  #-R ${REF}

done

