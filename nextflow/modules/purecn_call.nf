nextflow.enable.types = true

process PURECN_CALL {
    input:
    record(
        id: String,
        coverage: Path,
        mutect2_vcf: Path
    )
    normaldb: Path
    bait_intervals: String
    mapping_bias: Path
    snp_blacklist: Path
    default_purity: String
    default_ploidy: String

    output:
    record(
        id: id,
        call: file(output),
        purity: eval("head purity.txt"),
        ploidy: eval("head ploidy.txt"),
    )

    topic:
    record(id: id, process: "purecn_call", file: file("purecn_call.log")) >> "logs"

    script:
    //
    output = "${id}-Call"
    blacklist_flag = snp_blacklist ? " --snp-blacklist ${snp_blacklist}" : ""
    mb_flag = mapping_bias ? " --mapping-bias-file ${mapping_bias}" : ""
    get_purity_ploidy = """
    if [[ -e ${output}/${id}.csv ]]; then
         cut -d, -f 2 ${output}/${id}.csv | tail -n 1 > purity.txt
         cut -d, -f 3 ${output}/${id}.csv | tail -n 1 > ploidy.txt
    else
        echo ${default_purity} > purity.txt
        echo ${default_ploidy} > ploidy.txt
    fi
    """
    args = task.ext.args ?: ""
    """
    Rscript \$PURECN/PureCN.R \\
        ${args} \\
        ${mb_flag} \\
        ${blacklist_flag} \\
        --out . \\
        --sampleid ${id} \\
        --tumor "${coverage}" \\
        --vcf "${mutect2_vcf}" \\
        --normaldb "${normaldb}" \\
        --intervals "${bait_intervals}" \\
        --genome ${params.genome_build}

    mkdir _tmp
    mv ${id}* _tmp
    mv _tmp ${output}

    ${get_purity_ploidy}

    get_nextflow_log.bash purecn_call.log
    """
}
