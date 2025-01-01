process BCFTOOLS_STATS {
    ext version: params.bcftools_version

    publishDir "${meta.out}", mode:"copy", saveAs: params.saveFn
    publishDir "${meta.log}", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(vcf)
    val(target_intervals)
    val(module_number)

    output:
    tuple val(meta), path(output), emit: plots
    path("${output}/${output}.py"), emit: py
    path("*.log")
    //

    script:
    output = Utl.getName(module_number, meta, "Bcftools_stats")
    check = file("${meta.out}/${output}")
    region_flag = target_intervals != "" ? "--regions-file ${target_intervals}" : ""
    args = task.ext.args.join(" ")
    if (check.exists()) {
        """
        cp -r ${check} .
        ln -sr ${meta.log}/bcftools_stats.log .
        """
    } else {
        """
        bcftools index ${vcf}
        bcftools stats --split-by-ID \\
            ${region_flag} \\
            ${vcf} > data.vchk
        plot-vcfstats data.vchk -p ${output} --no-PDF
        mv ${output}/plot.py ${output}/${output}.py

        get_nextflow_log.bash bcftools_stats.log
        """
    }
    //
}
