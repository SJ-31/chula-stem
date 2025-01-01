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
    path("*.log")
    //

    script:
    output = Utl.getName(module_number, meta, "SigProfilerAssignment")
    exome_flag = is_exome ? " --exome " : ""
    exclude_flag = exclude_file != "" ? " --exclude_file ${exclude_file} " : ""
    args = task.ext.args.join(" ")
    check = file("$meta.out/$output")
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
            cp \$v vcfs
            gunzip vcfs/\$v
        done

        SigProfilerAssignmentWrapper.py -i vcfs \\
            ${args} \\
            -o ${output} \\
            ${exome_flag} \\
            ${exclude_flag}

        mv ${output}/Assignment_Solution/* ${output}
        rmdir ${output}/Assignment_Solution

        get_nextflow_log.bash sigprofilerassignment.log
        """
    }
}
