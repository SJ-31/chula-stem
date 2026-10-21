process VEP {
    ext version: "113"

    input:
    tuple val(meta), path(vcf)
    val reference
    val panel_of_normals

    output:

    //
    tuple val(meta), path(output), emit: vcf
    tuple val(meta), path ("*.html"), emit: report
    tuple val(meta), path(tsv), emit: tsv
    tuple val(meta), path("*.log"), emit: log

    script:
    //
    output = Utl.getName(meta, "VEP", "vcf.gz")
    tsv = Utl.getName(meta, "VEP", "tsv")
    logfile = Utl.getName(meta, "VEP", "log")
    html = Utl.getName(meta, "VEP_summary", "html")
    variant_class = meta.variant_class ? meta.variant_class : "small"
    // One of "sv" or "small"
    normal_flag = !params.tumor_only ? "-n ${meta.RGSM_normal}" : ""
    input = panel_of_normals ? "pon_filtered.vcf.gz" : vcf
    if (panel_of_normals) {
        pon_filter_command = """
        bcftools view -T ^${panel_of_normals} ${vcf} -O z -W -o pon_filtered.vcf.gz
        """
    }
    else {
        pon_filter_command = ""
    }
    args = task.ext.args.join(" ")
    """
    ${pon_filter_command}
    
    bcftools view ${input} -H | head > check.txt

    if [[ -s check.txt ]]; then

            vep --cache \\
                ${args} \\
                --input_file ${input} \\
                --dir_cache ${task.ext.cache} \\
                --species ${task.ext.species} \\
                --fasta ${reference} \\
                --stats_file ${html} \\
                --hgvs \\
                --hgvsg \\
                --vcf \\
                --vcf_info_field ANN \\
                --compress_output bgzip \\
                --output_file ${output}

            format_vep_vcf -i ${output} -o ${tsv} -t ${params.source_name} \\
                -c ${variant_class} -v ANN -r ${meta.RGSM_tumor} ${normal_flag}

            get_nextflow_log.bash ${logfile}
    else
        touch ${output}
        touch ${html}
        echo '${input} has no variants!' > ${logfile}
        touch ${tsv}
    fi
    """
    // use --synonyms flag to unify the naming systems (can just use your existing renaming files)
}
