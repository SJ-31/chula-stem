process DELLY_SV {
    ext version: "1.3.1"

    input:
    tuple val(meta), path(normal), path(tumor), path(indices, arity: "2")
    val reference
    val exclude

    output:

    //
    tuple val(meta), path(out), emit: variants
    tuple val(meta), path("*.log"), emit: logfile

    script:
    //
    out = Utl.getName(meta, "DellySV", "vcf.gz")
    args = task.ext.args.join(" ")
    exclude_flag = exclude ? "-x ${exclude}" : ""
    """
    get_sample_file () {
        normal_sample=${meta.RGSM_normal}
        tumor_sample=\$(bcftools query -l "\$1" | sed "s/\$normal_sample//" | xargs)
        echo -e "\${normal_sample}\tcontrol" > samples.tsv
        echo -e "\${tumor_sample}\ttumor" >> samples.tsv
    }

    delly call \
        ${args} \
        -g ${reference} \
        ${exclude_flag} \
        -o tmp.vcf \
        ${tumor} \
        ${normal}

    get_sample_file tmp.vcf

    delly filter -f somatic -o pre.bcf -s samples.tsv tmp.vcf # Attempts to filter false
    # positives and germline SVs

    # Genotype the somatic sites
    delly call ${args} \
        -g ${reference} \
        -v pre.bcf \
        ${exclude_flag} \
        -o geno.bcf \
        ${tumor} \
        ${normal}

    # Post-filter
    delly filter -f somatic -o tmp2.bcf -s samples.tsv geno.bcf

    bcftools view -s "${meta.RGSM_normal},${meta.RGSM_tumor}" tmp2.bcf | \
        vcf_info_add_tag.bash -n ${params.source_name} \
            -d "${params.source_description}" \
            -b '.' \
            -t String \
            -a delly \
            -o "${out}"

    get_nextflow_log.bash dellySV.log
    """
    //
}
