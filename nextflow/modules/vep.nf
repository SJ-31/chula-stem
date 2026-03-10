process VEP {
    ext version: "113"

    publishDir "$meta.out", mode: "copy", saveAs: params.saveFn
    publishDir "$meta.log", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(vcf)
    val(reference)
    val(panel_of_normals)
    val(module_number)
    //

    output:
    tuple val(meta), path(output), emit: vcf
    path("*.html"), emit: report
    tuple val(meta), path(tsv), emit: tsv
    path("*.log")
    //

    script:
    output = Utl.getName(module_number, meta, "VEP", "vcf.gz")
    tsv = Utl.getName(module_number, meta, "VEP", "tsv")
    logfile = Utl.getName(module_number, meta, "VEP", "log")
    html = Utl.getName(module_number, meta, "VEP_summary", "html")
    check = file("$meta.out/$output")
    check2 = file("$meta.out/$tsv")
    variant_class = meta.variant_class ? meta.variant_class : "small"
    // One of "sv" or "small"
    normal_flag = !params.tumor_only ? "-n ${meta.RGSM_normal}" : ""
    input = panel_of_normals ? "pon_filtered.vcf" : vcf
    if (panel_of_normals) {
        pon_filter_command = """
        rtg vcffilter -i ${vcf} --exclude-vcf=${panel_of_normals} -o pon_filtered.vcf
        """
    } else {
        pon_filter_command = ""
    }
    args = task.ext.args.join(" ")
    if (check.exists() && check2.exists()) {
        """
        ln -sr $check .
        ln -sr $check2 .
        ln -sr $meta.log/${logfile} .
        ln -sr $meta.out/${html} .
        """
    } else {
        """
        bcftools view ${input} -H | head > check.txt

        if [[ -s check.txt ]]; then
                ${pon_filter_command}

                vep --cache \\
                    $args \\
                    --input_file ${input} \\
                    --dir_cache ${task.ext.cache} \\
                    --species ${task.ext.species} \\
                    --fasta $reference \\
                    --stats_file ${html} \\
                    --hgvs \\
                    --hgvsg \\
                    --vcf \\
                    --vcf_info_field ANN \\
                    --compress_output bgzip \\
                    --output_file $output

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
    }
    // use --synonyms flag to unify the naming systems (can just use your existing renaming files)
}
