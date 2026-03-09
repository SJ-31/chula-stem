process PURECN_CALL {
    ext version: "2.16.1"

    publishDir "${meta.out}", mode:"copy", saveAs: params.saveFn
    publishDir "${meta.log}", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(coverage), path(mutect2_vcf)
    path(normaldb)
    path(bait_intervals)
    val(snp_blacklist)
    val(default_purity)
    val(default_ploidy)
    val(module_number)
    //

    output:
    tuple val(meta), path(output), emit: call
    tuple val(meta), eval("head purity.txt"), eval("head ploidy.txt"), emit: purity_ploidy
    path("*.log")
    //

    script:
    output = Utl.getName(module_number, meta, "Call", "")
    check = file("${meta.out}/${output}")
    blacklist_flag = snp_blacklist ? " --snp-blacklist ${snp_blacklist}" : ""
    sample_id = "${module_number}_${meta.id}"
    args = task.ext.args.join(" ")
    if (check.exists()) {
        """
        ln -sr ${check} .
        ln -sr ${meta.log}/purecn_coverage.log .
        """
    } else {
        """
        Rscript \$PURECN/PureCN.R \\
            ${args} \\
            ${blacklist_flag} \\
            --out . \\
            --sampleid ${sample_id} \\
            --tumor "${coverage}" \\
            --vcf "${mutect2_vcf}" \\
            --normaldb "${normaldb}" \\
            --intervals "${bait_intervals}" \\
            --genome ${params.genome_build}

        mkdir ${output}
        mv "${sample_id}*" ${output}

        if [[ -e ${output}/${sample_id}.csv ]]; then
            cut -d, -f 2 ${output}/${sample_id}.csv | tail -n 1 > purity.txt
            cut -d, -f 3 ${output}/${sample_id}.csv | tail -n 1 > ploidy.txt
        else
            echo ${default_purity} > purity.txt
            echo ${default_ploidy} > ploidy.txt
        fi

        get_nextflow_log.bash purecn_call.log
        """
    }
    //
}
