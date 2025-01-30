process DUPRADAR {
    ext version: params.bcftools_version

    publishDir "${meta.out}", mode:"copy", saveAs: params.saveFn
    publishDir "${meta.log}", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(bam)
    val(reference) // GTF reference
    val(strandedness)
    val(paired)
    val(module_number)

    output:
    tuple val(meta), path(output), emit: plots
    tuple val(meta), path("${output}/${output}.tsv"), emit: tsv
    path("*.log")
    //

    script:
    output = Utl.getName(module_number, meta, "Dupradar")
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
