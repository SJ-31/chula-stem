process PURECN_CALL {
    ext version: "2.16.1"

    publishDir "${meta.out}", mode:"copy", saveAs: params.saveFn
    publishDir "${meta.log}", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(coverage), path(mutect2_vcf)
    val(normaldb)
    val(bait_intervals)
    val(mapping_bias)
    val(snp_blacklist)
    val(default_purity)
    val(default_ploidy)
    val(module_number)
    //

    output:
    tuple val(meta), path(output), emit: call
    tuple val(meta), eval("head purity.txt"), eval("head ploidy.txt"), emit: purity_ploidy
    path("*.log"), optional: true
    //

    script:
    output = Utl.getName(module_number, meta, "Call", "")
    check = file("${meta.out}/${output}")
    blacklist_flag = snp_blacklist ? " --snp-blacklist ${snp_blacklist}" : ""
    mb_flag = mapping_bias ? " --mapping-bias-file ${mapping_bias}" : ""
    get_purity_ploidy = """
    if [[ -e ${output}/${meta.id}.csv ]]; then
         cut -d, -f 2 ${output}/${meta.id}.csv | tail -n 1 > purity.txt
         cut -d, -f 3 ${output}/${meta.id}.csv | tail -n 1 > ploidy.txt
    else
        echo ${default_purity} > purity.txt
        echo ${default_ploidy} > ploidy.txt
    fi
    """
    args = task.ext.args.join(" ")
    if (check.exists()) {
        """
        cp -r ${check} .
        ln -sr ${meta.log}/purecn_call.log .

        ${get_purity_ploidy}
        """
    } else {
        """
        Rscript \$PURECN/PureCN.R \\
            ${args} \\
            ${mb_flag} \\
            ${blacklist_flag} \\
            --out . \\
            --sampleid ${meta.id} \\
            --tumor "${coverage}" \\
            --vcf "${mutect2_vcf}" \\
            --normaldb "${normaldb}" \\
            --intervals "${bait_intervals}" \\
            --genome ${params.genome_build}

        mkdir _tmp
        mv ${meta.id}* _tmp
        mv _tmp ${output}

        ${get_purity_ploidy}

        get_nextflow_log.bash purecn_call.log
        """
    }
    //
}
