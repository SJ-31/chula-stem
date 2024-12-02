process VEP {
    ext version: "113"

    publishDir "$meta.out", mode: "copy", saveAs: params.saveFn
    publishDir "$meta.log", mode: "copy", pattern: "*.{log,html}"

    input:
    tuple val(meta), path(vcf)
    val(reference)
    val(module_number)
    //

    output:
    tuple val(meta), path(output), emit: vcf
    tuple val(meta), path("*.html"), emit: report
    path("*.log")
    //

    script:
    output = Utils.getName(module_number, meta, "VEP", "vcf.gz")
    tsv = Utils.getName(module_number, meta, "VEP", "tsv")
    check = file("$meta.out/$output")
    args = task.ext.args.join(" ")
    if (check.exists()) {
        """
        ln -sr $check .
        ln -sr $meta.log/vep.log .
        ln -sr $meta.log/vep_stats.html .
        """
    } else {
        """
        vep --cache \\
            $args \\
            --input_file $vcf \\
            --dir_cache ${task.ext.cache} \\
            --species ${task.ext.species} \\
            --fasta $reference \\
            --stats_file vep_stats.html \\
            --vcf \\
            --vcf_info_field ANN \\
            --compress_output bgzip \\
            --output_file $output

        format_vep_vcf -i ${output} -o ${tsv} -t ${params.source_name} -v ANN

        get_nextflow_log.bash vep.log
        """
    }
    // use --synonyms flag to unify the naming systems (can just use your existing renaming files)
}
