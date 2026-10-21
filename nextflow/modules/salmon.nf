process SALMON {
    ext version: "1.10.3"

    input:
    tuple val(meta), val(reads_or_bam)
    val index // Salmon index
    val reference // Must be a reference transcriptome that reads were aligned to
    // ignored if using reads mode
    val gene_map // GTF/GFF file mapping transcripts to genes
    val strandedness // foward, reverse or null
    val relative_orientation // inward,outward or matching

    output:

    //
    tuple val(meta), path(sf), emit: transcripts
    tuple val(meta), path(sf_genes), emit: genes
    tuple val(meta), path("*.log"), emit: log

    script:
    //
    sf = Utl.getName(meta, "salmon", "sf")
    sf_genes = Utl.getName(meta, "salmon.genes", "sf")
    logfile = Utl.getName(meta, "salmon", "log")
    args = task.ext.args.join(" ")

    if (relative_orientation == "inward") {
        libtype = "I"
    }
    else if (relative_orientation == "outward") {
        libtype = "O"
    }
    else if (relative_orientation == "matching") {
        libtype = "M"
    }
    else {
        throw new Exception("relative orientation not recognized!")
    }

    if (strandedness == "forward") {
        libtype = "${libtype}SF"
    }
    else if (strandedness == "reverse") {
        libtype = "${libtype}SR"
    }
    else {
        libtype = "${libtype}U"
    }
    if (index) {
        """
        salmon quant ${args} \\
            -i ${index} \\
            -l ${libtype} \\
            -1 ${reads_or_bam[0]} -2 ${reads_or_bam[1]} \\
            -g ${gene_map} \\
            -o salmon_quant

        mv salmon_quant/quant.sf ${sf}
        mv salmon_quant/quant.genes.sf ${sf}
        get_nextflow_log.bash ${logfile}
        """
    }
    else {
        """
        salmon quant -l ${libtype} \\
            ${args} \\
            -a ${reads_or_bam} \\
            -t ${reference}  \\
            -g ${gene_map} \\
            -o salmon_quant

        mv salmon_quant/quant.sf ${sf}
        mv salmon_quant/quant.genes.sf ${sf}
        get_nextflow_log.bash ${logfile}
        """
    }
    //
}
