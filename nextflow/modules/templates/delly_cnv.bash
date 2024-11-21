get_sample_file () {
    normal_sample=!{meta.RGSM_normal}
    tumor_sample=$(bcftools query -l "$1" | sed "s/$normal_sample//" | xargs)
    echo -e "${normal_sample}\tcontrol" > samples.tsv
    echo -e "${tumor_sample}\ttumor" >> samples.tsv
}

# Segment tumor genome
# TODO: can't fill these parameters until you have software that estimates ploidy and purity
# purity=$()
# ploidy=
delly cnv --segmentation \
    --cnv-size "TODO" \
    --cn_offset "TODO" \
    --sdrd "TODO" \
    --purity "TODO" \
    --outfile tumor.bcf \
    --covfile !{covfile} \
    --genome !{reference} \
    --mappability !{mappability} \
    !{tumor}
# Detection sensitivity is controlled by
# --cnv-size (-z):
# --cn_offset (-t):
#   Is the minimum copy-number shift for segmentation and somatic classification
#   Can calculate as purity * ploidy + 2 (diploid number) * (1-purity)
# --sdrd (read-depth shift, -x):
#
# which should be determined dynamically for each tumor sample

# Genotype somatic CNVs in normal sample
delly cnv --segmentation \
    --vcffile tumor.bcf \
    --outfile normal.bcf \
    --genome !{reference} \
    --mappability \
    !{normal}

bcftools merge -m id -O b -o tmp.bcf tumor.bcf normal.bcf
get_sample_file tmp.bcf

delly classify --pass --filter somatic --outfile !{out} -s samples.tsv tmp.bcf

get_nextflow_log.bash dellyCNV.log
