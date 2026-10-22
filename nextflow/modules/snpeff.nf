process SNPEFF {
    ext version: "5.2e"
    
    
    

    input:
    tuple val(meta), path(vcf)
    val(reference)
    val(is_somatic)
    
    //

    output:
    tuple val(meta), path("${output}.gz"), emit: vcf
    tuple val(meta), path("*.html"), path("*.txt"), emit: report
    tuple val(meta), path("*.log"), emit: logfile
    //

    script:
    output = Utl.getName(meta, "snpeFF", "vcf")
    report = Utl.getName(meta, "snpEff_summary", "html")
    genes_file = Utl.getName(meta, "snpEff_genes", "txt")
    args = task.ext.args.join(" ")
        """
        if [[ ${is_somatic} == "true" ]]; then
            # This only works if the given vcf has only one patient
            normal_sample=${meta.RGSM_normal}
            cancer_sample=$(bcftools query -l ${vcf} | sed "s/\$normal_sample//" | xargs)
            echo -e "\${normal_sample}\t\${cancer_sample}" > cancer_samples.txt

            ${params.snpEff} ${args} \
                -nodownload  \
                -cancerSamples cancer_samples.txt \
                ${reference} \
                ${vcf} > ${output}
        else
            ${params.snpEff} ${args} \
                -nodownload  \
                ${reference} \
                ${vcf} > ${output}
        fi


        bgzip ${output}

        mv snpEff_summary.html ${report}
        mv snpEff_genes.txt ${genes_file}
        get_nextflow_log.bash snpEff.log
        """
}
