process SALMON {
    ext version: "1.10.3"

    publishDir "${meta.out}", mode:"copy", saveAs: params.saveFn
    publishDir "${meta.log}", mode: "copy", pattern: "*.{log,json}"

    input:
    tuple val(meta), val(reads_or_bam)
    val(index) // Salmon index
    val(reference) // Must be a reference transcriptome that reads were aligned to
    val(gene_map) // GTF/GFF file mapping transcripts to genes
    val(strandedness) // foward, reverse or null
    val(relative_orientation) // inward,outward or matching
    val(module_number)
    //

    output:
    tuple val(meta), path(sf), emit: transcripts
    tuple val(meta), path(sf_genes), emit: genes
    path("*.log")
    //

    script:
    sf = Utl.getName(module_number, meta, "salmon", "sf")
    sf_genes = Utl.getName(module_number, meta, "salmon.genes", "sf")
    logfile = Utl.getName(module_number, meta, "salmon", "log")
    checks = [sf, sf_genes].collect({ file("${meta.out}/${it}") })
    args = task.ext.args.join(" ")

    if (relative_orientation == "inward") {
        libtype = "I"
    } else if (relative_orientation == "outward") {
        libtype = "O"
    } else if (relative_orientation == "matching") {
        libtype = "M"
    } else {
        throw new Exception("relative orientation not recognized!")
    }

    if (strandedness == "forward") {
        libtype = "${libtype}SF"
    } else if (strandedness == "reverse") {
        libtype = "${libtype}SR"
    } else {
        libtype = "${libtype}U"
    }
    if (checks.every({ it.exists() })) {
        """
        ln -sr ${checks[0]} .
        ln -sr ${checks[1]} .
        ln -sr ${meta.log}/${logfile} .
        """
    } else if (index) {
        """
        salmon quant ${args} \\
            -i ${index} \\
            -l ${libtype} \\
            -1 ${reads_or_bam[0]} -2 ${reads_or_bam[1]} \\
            -g ${gene_map}

        mv salmon_quant/quant.sf ${sf}
        mv salmon_quant/quant.genes.sf ${sf}
        get_nextflow_log.bash ${logfile}
        """
    } else {
        """
        salmon quant -l ${libtype} \\
            ${args} \\
            -a ${reads_or_bam} \\
            -t ${reference}  \\
            -g ${gene_map}

        mv salmon_quant/quant.sf ${sf}
        mv salmon_quant/quant.genes.sf ${sf}
        get_nextflow_log.bash ${logfile}
        """
    }
    //
}
