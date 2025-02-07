process COMBINE_COUNTS {
    ext version: "1"

    publishDir "${meta.out}", mode:"copy", saveAs: params.saveFn
    publishDir "${meta.log}", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), val(file_data)
    val(module_number)
    //

    output:
    tuple val(meta), path(output)
    path("*.log")
    //

    script:
    output = Utl.getName(module_number, meta, "All_counts", "tsv")
    check = file("${meta.out}/${output}")
    csv_rep = Utl.listOfMaps2Csv(file_data)
    if (check.exists()) {
        """
        ln -sr ${check} .
        ln -sr ${meta.log}/combine_counts.log .
        """
    } else {
        """
        echo -e "${csv_rep}" > spec.csv

        cli.R -c combine_counts \\
            -i spec.csv -o ${output}

        get_nextflow_log.bash combine_counts.log
        """
    }
    //
}
