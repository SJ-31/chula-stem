process GRIDSS {
    ext version: "2.13.2"

    publishDir "$meta.out", mode:"copy", saveAs: params.saveFn
    publishDir "$meta.log", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(normal), path(tumor), path(indices, arity: "2")
    val(reference)
    val(blacklist) // list of genomic regions to exclude, gridss docs recommends using
    // ENCODE blacklists
    val(module_number)
    //

    output:
    tuple val(meta), path(filtered), emit: variants
    path(all)
    path("*.log")
    //

    shell:
    prefix = "${module_number}-${meta.filename}"
    all = "${prefix}-Gridss_all.vcf.gz"
    filtered = "${prefix}-Gridss_confident.vcf.gz"
    check1 = file("${meta.out}/${all}")
    check2 = file("${meta.out}/${filtered}")
    args = task.ext.args.join(" ")
    if (check1.exists() && check2.exists()) {
        '''
        ln -sr !{check2} .
        ln -sr !{check1} .
        ln -sr !{meta.log}/gridss.log .
        '''
    } else {
        '''
        gridss \\
            --threads !{task.ext.threads} \\
            --reference !{reference} \\
            --output tmp.vcf \\
            !{normal} \\
            !{tumor}

        gridss_somatic_filter \\
            --input tmp.vcf \\
            --output confident.vcf \\
            --fulloutput all.vcf

        names=("!{filtered}" "!{all}")
        files=(confident.vcf.bgz all.vcf.bgz)

        for i in $(seq 0 1); do
            name="${names[i]}"
            file="${files[i]}"

            if [[ "${file}" == "confident.vcf.bgz" ]]; then
                bcftools view -s "!{meta.RGSM_normal},!{meta.RGSM_tumor}" -O z "${file}" > tmp.vcf.gz
                f="tmp.vcf.gz"
            else
                f="${file}"
            fi

            vcf_info_add_tag.bash.r -n !{params.source_name} \\
                -d '!{params.source_description}' \\
                -b '.' \\
                -t String \\
                -a gridss \\
                -i "${f}" \\
                -o "$name"
        done

        get_nextflow_log.bash gridss.log
        '''
    }
    //
    // TOOD: you can supply a pon with this
}
