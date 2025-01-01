process DELLY_CNV {
    ext version: "1.3.1"
 
    publishDir "$meta.out", mode: "copy", saveAs: params.saveFn
    publishDir "$meta.log", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(normal), path(tumor), path(indices, arity: "2"), path(covfile), val(purity), val(ploidy)
    val(reference)
    val(mappability)
    val(module_number)
    //

    output:
    tuple val(with_caller), path(out)
    path(segmentation)
    path("*.log")
    //

    script:
    out = Utl.getName(module_number, meta, "DellyCNV", ".vcf.gz")
    segmentation = Utl.getName(module_number, meta, "DellySegmentation", "vcf.gz")
    with_caller = meta + ["caller": "dellyCNV"]
    check = file("${meta.out}/${out}")
    check2 = file("${meta.out}/${segmentation}")
    args = task.ext.args.join(" ")
    if (check.exists() && check2.exists()) {
        """
        ln -sr ${check} .
        ln -sr ${check2} .
        ln -sr ${meta.log}/dellyCNV.log .
        """
    } else {
        """
        get_sample_file () {
            normal_sample=${meta.RGSM_normal}
            tumor_sample=\$(bcftools query -l "\$1" | sed "s/\$normal_sample//" | xargs)
            echo -e "\${normal_sample}\tcontrol" > samples.tsv
            echo -e "\${tumor_sample}\ttumor" >> samples.tsv
        }

        # Segment tumor genome
        delly cnv --segmentation \
            ${args} \
            --ploidy ${ploidy} \
            --purity ${purity} \
            --outfile tumor.bcf \
            --covfile ${covfile} \
            --genome ${reference} \
            --mappability ${mappability} \
            ${tumor}
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
            --genome ${reference} \
            --mappability ${mappability} \
            ${normal}

        bcftools merge -m id -O b -o tmp.bcf tumor.bcf normal.bcf
        bcftools index tmp.bcf
        get_sample_file tmp.bcf

        delly classify --pass --filter somatic --outfile tmp.vcf.gz -s samples.tsv tmp.bcf

        bcftools view -s "${meta.RGSM_normal},${meta.RGSM_tumor}" tmp.vcf.gz | \
            vcf_info_add_tag.bash -n ${params.source_name} \
                -d "${params.source_description}" \
                -b '.' \
                -t String \
                -a dellyCNV \
                -o ${out}

        mv tumor.bcf "${segmentation}"

        get_nextflow_log.bash dellyCNV.log
        """
    }
    //
}
