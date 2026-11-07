nextflow.enable.types = true

process SIGPROFILERASSIGNMENT {
    
    input:
    record(id: String,
           vcfs: Bag<Path>)
    is_exome: Boolean
    exclude_file: Path

    output:
    record(id: id,
           dir: file(all),
           activities: activities_file
    )

    topic:
    record(id: id, process: "sigprofilerassignment_collect", file: file("sigprofilerassignment.log")) >> "logs"

    script:
    //
    output = "${id}-SigProfilerAssignment"
    exome_flag = is_exome ? " --exome " : ""
    exclude_flag = exclude_file != "" ? " --exclude_file ${exclude_file} " : ""
    args = task.ext.args ?: ""
    activities_file = "${output}/Activities/Assignment_Solution_Activities.txt"
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

    if [[ ! -z \$(ls -A ${output}/Assignment_Solution) ]]; then
        mv ${output}/Assignment_Solution/* ${output}
    fi
    
    rmdir ${output}/Assignment_Solution

    if [[ ! -e ${output}/Activities ]]; then
        mkdir ${output}/Activities
    fi

    if [[ ! -e '${activities_file}' ]]; then
        echo EMPTY > ${activities_file}
    fi

    get_nextflow_log.bash sigprofilerassignment.log
    """
}

process SIGPROFILERASSIGNMENT_COLLECT {
    input:
    tuple val(meta), path(activities, stageAs: "?/*")

    output:

    //
    path output

    script:
    //
    output = Utl.getName(meta, "SigProfilerAssignment_activities", "tsv")
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
    //
}
