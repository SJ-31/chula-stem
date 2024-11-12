process VEP {
    ext version: "113"

    publishDir "$meta.out", mode: "copy"
    publishDir "$meta.log", mode: "copy", pattern: "*.{log,html}"

    input:
    tuple val(meta), path(vcf)
    val(reference)
    val(module_number)
    //

    output:
    tuple val(meta), path(output)
    path("*.log")
    path("*.html")
    //

    script:
    output = "${module_number}-${meta.id}_vep.vcf.gz"
    check = file("$meta.out/$output")
    def args = task.ext.args.join(" ")
    if (check.exists()) {
        """
        cp $check .
        cp $meta.log/vep.log .
        cp $meta.log/vep_stats.html .
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
            --compress_output \\
            --output_file $output

        cp .command.out vep.log
        """
    }
    //
}
