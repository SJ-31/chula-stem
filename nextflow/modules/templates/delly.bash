get_sample_file () {
    normal_sample=!{meta.RGSM_normal}
    tumor_sample=$(bcftools query -l "$1" | sed "s/$normal_sample//" | xargs)
    echo -e "${normal_sample}\tcontrol" > samples.tsv
    echo -e "${tumor_sample}\ttumor" >> samples.tsv
}

delly call \
    !{args} \
    -g !{reference} \
    !{exclude_flag} \
    -o tmp.vcf \
    !{tumor} \
    !{normal}

get_sample_file tmp.vcf

delly filter -f somatic -o pre.bcf -s samples.tsv tmp.vcf # Attempts to filter false
# positives and germline SVs

# Genotype the somatic sites
delly call !{args} \
    -g !{reference} \
    -v pre.bcf \
    !{exclude_flag} \
    -o geno.bcf \
    !{tumor} \
    !{normal}

# Post-filter
delly filter -f somatic -o !{out} -s samples.tsv geno.bcf

get_nextflow_log.bash delly.log
