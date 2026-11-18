nextflow.enable.types = true

// WARNING: [2026-11-18 Wed] sigprofilerassignment needs to be installed
// with matplotlib<=3.9
// e.g. conda create -n SigProfiler sigprofilerassignment matplotlib=3.9 -c bioconda
process SIGPROFILERASSIGNMENT {
    input:
    record(
        id: String,
        vcfs: List<Path>
    )
    is_exome: Boolean
    exclude_file: Path?
    volume: String // path to the installed sigprofiler reference
    // genome data

    output:
    record(
        id: id,
        dir: file(output),
        activities: file(activities_file),
    )

    topic:
    record(id: id, process: "sigprofilerassignment_collect", file: file("sigprofilerassignment.log")) >> "logs"

    script:
    //
    output = "${id}-SigProfilerAssignment"
    exome_flag = is_exome ? " --exome " : ""
    exclude_flag = exclude_file ? " --exclude_file ${exclude_file} " : ""
    args = task.ext.args ?: ""
    activities_file = "${output}/Activities/Assignment_Solution_Activities.txt"
    """
    mkdir vcfs

    for v in *vcf.gz; do
        bcftools index \$v
        sigprofiler_assignment_cleanup.bash \$v vcfs
    done

    SigProfilerAssignmentWrapper.py -i vcfs \\
        ${args} \\
        -o ${output} \\
        -l ${volume} \\
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
    activities: Bag<Path>
    sample_pattern: String

    output:
    file(out)

    script:
    //
    out = "SigProfilerAssignment_activities.tsv"
    """
    #!/usr/bin/env Rscript

    library(tidyverse)
    files <- list.files(pattern = "Assignment_Solution_Activities.txt",
        recursive = TRUE)
    combined <- lapply(files, read_tsv) |> bind_rows()
    if (nrow(combined) > 0) {
        combined[['Samples']] <- map_chr(combined[['Samples']], function(x) {
    str_extract(x, "${sample_pattern}", group = 1)
        })
    }
    write_tsv(combined, "${out}")
    """
    // TODO: you need to modify the sample pattern depending on the
    // filename
}
