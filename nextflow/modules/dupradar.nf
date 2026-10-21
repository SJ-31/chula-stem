process DUPRADAR {
    ext version: params.bcftools_version

    
    

    input:
    tuple val(meta), path(bam)
    val(reference) // GTF reference
    val(strandedness)
    val(paired)
    

    output:
    tuple val(meta), path(output), emit: plots
    tuple val(meta), path("${output}/${output}.tsv"), emit: tsv
    tuple val(meta), path("*.log"), emit: log
    //

    script:
    output = Utl.getName(meta, "Dupradar")
    check = file("${meta.out}/${output}")
    paired_flag = paired ? "" : " --unpaired "
    if (check.exists()) {
        """
        cp -r ${check} .
        ln -sr ${meta.log}/dupradar.log .
        """
    } else {
        """
        mkdir ${output}
        dupradar.R -b ${bam} \\
            -g ${reference} \\
            -o ${output}/${output}.tsv \\
            -s ${strandedness} \\
            --box_plot ${output}/box_plot.png \\
            --density_plot ${output}/density_plot.png \\
            ${paired_flag}

        get_nextflow_log.bash dupradar.log
        """
    }
    //
}
