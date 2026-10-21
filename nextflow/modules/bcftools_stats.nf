process BCFTOOLS_STATS {
    ext version: params.bcftools_version

    
    

    input:
    tuple val(meta), path(vcf)
    val(target_intervals)
    

    output:
    tuple val(meta), path(output), emit: plots
    path("${output}/${output}.py"), emit: py
    tuple val(meta), path("*.log"), emit: log
    //

    script:
    output = Utl.getName(meta, "Bcftools_stats")
    check = file("${meta.out}/${output}")
    region_flag = target_intervals != "" ? "--regions-file ${target_intervals}" : ""
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
