process SIGPROFILERASSIGNMENT {
    ext version: "0.1.9"

    publishDir "$meta.out", mode: "copy", saveAs: params.saveFn
    publishDir "$meta.log", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(vcfs)
    val(is_exome)
    val(exclude_file)
    val(module_number)
    //

    output:
    tuple val(meta), path(output), emit: data
    tuple val(meta), path(activities_file), emit: activities
    path("*.log")
    //

    script:
    output = Utl.getName(module_number, meta, "SigProfilerAssignment")
    exome_flag = is_exome ? " --exome " : ""
    exclude_flag = exclude_file != "" ? " --exclude_file ${exclude_file} " : ""
    args = task.ext.args.join(" ")
    check = file("$meta.out/$output")
    activities_file = "${output}/Activities/Assignment_Solution_Activities.txt"
    if (check.exists()) {
        """
        cp -r $check .
        ln -sr $meta.log/sigprofilerassignment.log .
        """
    } else {
        """
        export SIGPROFILERASSIGNMENT_VOLUME=${task.ext.volume}

        mkdir vcfs

        for v in *vcf.gz; do
            sigprofiler_assignment_cleanup.bash \$v vcfs
        done

        SigProfilerAssignmentWrapper.py -i vcfs \\
            ${args} \\
            -o ${output} \\
            ${exome_flag} \\
            ${exclude_flag}

        mv ${output}/Assignment_Solution/* ${output}
        rmdir ${output}/Assignment_Solution

        if [[ ! -e '${activities_file}' ]]; then
            echo EMPTY > ${activities_file}
        fi

        get_nextflow_log.bash sigprofilerassignment.log
        """
    }
}

process SIGPROFILERASSIGNMENT_COLLECT {
    publishDir "${meta.out}", mode:"copy", saveAs: params.saveFn

    input:
    tuple val(meta), path(activities, stageAs: "?/*")
    val(module_number)
    //

    output:
    path(output)
    //

    script:
    output = Utl.getName(module_number, meta, "SigProfilerAssignment_activities", "tsv")
    check = file("${meta.out}/${output}")
    if (check.exists()) {
        """
        ln -sr ${check} .
        """
    } else {
        """
        #!/usr/bin/env Rscript

        library(tidyverse)
        files <- list.files(pattern = "Assignment_Solution_Activities.txt",
            recursive = TRUE)
        combined <- lapply(files, read_tsv) |> bind_rows()
        if (nrow(combined) > 0) {
            combined[['Samples']] <- map_chr(combined[['Samples']], function(x) {
            str_extract(x, '[0-9]+-(.*)-Small_high_conf', group = 1)
            })
        }
        write_tsv(combined, "${output}")
        """
    }
    //
}
