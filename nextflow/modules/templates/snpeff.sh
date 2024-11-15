if [[ !{is_somatic} == "true" ]]; then
    # This only works if the given vcf has only one patient
    normal_sample=!{meta.RGSM_normal}
    cancer_sample=$(bcftools query -l !{vcf} | sed "s/$normal_sample//" | xargs)
    echo -e "${normal_sample}\t${cancer_sample}" > cancer_samples.txt

    !{params.snpEff} !{args} \
        -nodownload  \
        -cancerSamples cancer_samples.txt \
        !{reference} \
        !{vcf} > !{output}
else
    !{params.snpEff} !{args} \
        -nodownload  \
        !{reference} \
        !{vcf} > !{output}
fi


bgzip !{output}

mv snpEff_summary.html !{report}
mv snpEff_genes.txt !{genes_file}
get_nextflow_log.bash snpEff.log
