process COMBINE_COUNTS {
    ext version: "1"

    
    

    input:
    tuple val(meta), val(file_data)
    
    //

    output:
    tuple val(meta), path(output)
    tuple val(meta), path("*.log"), emit: log
    //

    script:
    output = Utl.getName(meta, "All_counts", "tsv")
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
