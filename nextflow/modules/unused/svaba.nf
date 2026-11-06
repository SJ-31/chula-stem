process SVABA {
    ext version: "1.2.0"

    input:
    tuple val(meta), path(normal), path(tumor), path(indices, arity: "2")
    path reference
    path blacklist

    output:

    // [2026-10-22 Thu] TODO: specify the output correctly
    tuple val(meta), path("*_Manta.vcf.gz"), emit: variants
    tuple val(meta), path("*.log"), emit: logfile

    script:
    out = Utl.getName(meta, "SvABA")
    args = task.ext.args.join(" ")
    """
    svaba run \\
        -n ${normal} \\
        -t ${tumor} \\
        -G ${reference} \\
        -a ${out} \\
        --blacklist ${blacklist} \\
        ${args}

    svaba postprocess -i ${out} ${args}
    
    svaba tovcf -i ${out}.bps.sorted.dedup.txt.gz -b ${tumor} -a ${out}
    
    get_nextflow_log.bash svaba.log
    """
}
