process MERGE_VCFS {
    ext version: params.gatk_version

    
    

    input:
    tuple val(meta), path(vcfs)
    
    //

    output:
    tuple val(meta), path(output)
    tuple val(meta), path("*.log"), emit: log
    //

    script:
    output = Utl.getName(meta, "All_annotations", "vcf.gz")
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
