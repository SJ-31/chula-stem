get_sample_file () {
    normal_sample=!{meta.RGSM_normal}
    tumor_sample=$(bcftools query -l "$1" | sed "s/$normal_sample//" | xargs)
    echo -e "${normal_sample}\tcontrol" > samples.tsv
    echo -e "${tumor_sample}\ttumor" >> samples.tsv
}

# Segment tumor genome
delly cnv --segmentation \
    !{args} \
    --ploidy !{ploidy} \
    --purity !{purity} \
    --outfile tumor.bcf \
    --covfile !{covfile} \
    --genome !{reference} \
    --mappability !{mappability} \
    !{tumor}
# Detection sensitivity is controlled by
# --cnv-size (-z):
# --cn_offset (-t):
#   Is the minimum copy-number shift for segmentation and somatic classification
# --sdrd (read-depth shift, -x):
#
# which should be determined dynamically for each tumor sample

# Genotype somatic CNVs in normal sample
delly cnv --segmentation \
    --vcffile tumor.bcf \
    --outfile normal.bcf \
    --genome !{reference} \
    --mappability !{mappability} \
    !{normal}

bcftools merge -m id -O b -o tmp.bcf tumor.bcf normal.bcf
bcftools index tmp.bcf
get_sample_file tmp.bcf

delly classify --pass --filter somatic --outfile tmp.vcf.gz -s samples.tsv tmp.bcf

vcf_info_add_tag.bash -n !{params.source_name} \
    -d "!{params.source_description}" \
    -b '.' \
    -t String \
    -a dellyCNV \
    -i tmp.vcf.gz \
    -o tmp2.vcf

bcftools view -s "!{meta.RGSM_normal},!{meta.RGSM_tumor}" -O z tmp2.vcf > !{out}

mv tumor.bcf "!{segmentation}"
bgzip "!{segmentation}"

get_nextflow_log.bash dellyCNV.log
