nextflow.enable.types = true

// Combine CNV output from multiple callers
// Annotate concordant CNV calls with ClassifyCNV

process MERGE_CLASSIFY_CNV {
    input:
    record(
        id: String,
        caller2cnv_file: Map<String, Path>
        )
    genome_build: String

    output:
    // Input format for the BED file is to have the columns below
    // 1. chromosome
    // 2. CNV start
    // 3. CNV end
    // 4. CNV type (DEL or DUP)
    record(id: id, concordant: file(concordant), discordant: file(discordant), classified: file(tsv, optional: true))

    topic:
    record(id: id, process: "merge_classify_cnv", file: file("merge_classify_cnv.log")) >> "logs"

    script:
    //
    input_str = caller2cnv_file.entrySet().collect { e -> "-i ${e.key}:${e.value}" }.join(" ")
    tsv = "${id}-ClassifyCNV.tsv"
    concordant = "${id}-CNV_concordant.csv"
    discordant = "${id}-CNV_discordant.csv"
    intermediate = "${id}_intermediate"
    """
    merge_cnv_calls ${input_str} -c ${concordant} \\
        -d ${discordant}

    xan select chr,start,end,cn ${concordant} | \\
        xan filter 'cn !=2' | \\
        xan map 'if(cn > 2, "DUP", "DEL") as cat' | \\
        xan drop cn | \\
        xan fmt -t '\t' | \\
        xan behead > tmp.bed

    if [[ -s tmp.bed ]]; then
        ClassifyCNV.py --infile tmp.bed \\
            --GenomeBuild ${genome_build} \\
            --outdir ${intermediate}
        mv "${intermediate}/Scoresheet.txt" ${tsv}
    fi

    get_nextflow_log.bash merge_classify_cnv.log
    """
}
