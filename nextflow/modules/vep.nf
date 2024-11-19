process VEP {
    ext version: "113"

    publishDir "$meta.out", mode: "copy", saveAs: { x -> x ==~ /.*\.log/ ? null : x }
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
    output = "${module_number}-${meta.filename}-vep.vcf.gz"
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

        get_nextflow_log.bash vep.log
        """
    }
    //
}
