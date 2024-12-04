process VEP {
    ext version: "113"

    publishDir "$meta.out", mode: "copy", saveAs: params.saveFn
    publishDir "$meta.log", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(vcf)
    val(reference)
    val(module_number)
    //

    output:
    tuple val(meta), path(output), emit: vcf
    path("*.html"), emit: report
    path(tsv)
    path("*.log")
    //

    script:
    output = Utils.getName(module_number, meta, "VEP", "vcf.gz")
    tsv = Utils.getName(module_number, meta, "VEP", "tsv")
    html = Utils.getName(module_number, meta, "VEP_summary", "html")
    check = file("$meta.out/$output")
    check2 = file("$meta.out/$tsv")
    args = task.ext.args.join(" ")
    if (check.exists() && check2.exists()) {
        """
        ln -sr $check .
        ln -sr $check2 .
        ln -sr $meta.log/vep.log .
        ln -sr $meta.out/${html} .
        """
    } else {
        """
        vep --cache \\
            $args \\
            --input_file $vcf \\
            --dir_cache ${task.ext.cache} \\
            --species ${task.ext.species} \\
            --fasta $reference \\
            --stats_file ${html} \\
            --vcf \\
            --vcf_info_field ANN \\
            --compress_output bgzip \\
            --output_file $output

        format_vep_vcf -i ${output} -o ${tsv} -t ${params.source_name} \\
            -v ANN -r ${meta.RGSM_tumor}

        get_nextflow_log.bash vep.log
        """
    }
    // use --synonyms flag to unify the naming systems (can just use your existing renaming files)
}
