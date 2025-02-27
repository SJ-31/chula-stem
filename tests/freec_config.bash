make_freec_config --chr_len_file /data/project/stemcell/shannc/reference/genomes/GRCh38.p14_filtered.fna.fai \
    --read_type paired_end \
    --normal /data/project/stemcell/shannc/tests/wes_test2/patient_10_normal/4-patient_10_normal-recal.bam \
    --tumor /data/project/stemcell/shannc/tests/wes_test2/patient_10_cancer/4-patient_10_cancer-recal.bam \
    --chr_files /data/project/stemcell/shannc/tests/wes_test/controlfreec/chrs \
    --snps /data/project/stemcell/shannc/reference/variants/dbSNP_renamed.vcf.gz \
    --intervals /data/project/stemcell/shannc/reference/test_target_intervals_sorted.bed \
    --output sample_config \
    --noisy_data
