process MERGE_VCFS {
    ext version: "4.6.1.0"

    publishDir "$meta.out", mode: "copy"
    publishDir "$meta.log", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(vcfs)
    val(module_number)
    //

    output:
    tuple val(meta), path(output)
    path("*.log")
    //

    script:
    output = "${module_number}-${meta.id}_AllAnnotations.vcf.gz"
    check = file("${meta.out}/${output}")
    args = task.ext.args.join(" ")
    def file_list = []
    for (p in vcfs) {
        file_list.add("-I ${p.baseName}")
    }
    def merge_flag = file_list.join(" ")
    if (check.exists()) {
        """
        ln -sr $check .
        ln -sr ${meta.log}/vcf_merge_annotations.log .
        """
    } else {
        """
        gatk MergeVcfs ${merge_flag} \\
            -O $check

        get_nextflow_log.bash vcf_merge_annotations.log
        """
    }
    //
}
