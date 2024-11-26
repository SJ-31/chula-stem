process CONCAT_VCF {
    ext version: "1.21"
    // Concatenate groups vcf files of the same samples

    publishDir "$meta.out", mode: "copy", saveAs: params.saveFn
    publishDir "$meta.log", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(vcfs)
    val(suffix) // Suffix for ouptut file
    val(module_number)
    //

    output:
    tuple val(meta), path(output), emit: vcf
    path("*.log")
    //

    shell:
    output = "${module_number}-${meta.filename}-${suffix}.vcf.gz"
    check = file("${meta.out}/${output}")
    if (check.exists()) {
        '''
        ln -sr !{check} .
        ln -sr !{meta.log}/concat_vcf.log .
        '''
    } else {
        '''
        for i in *.vcf.gz; do
            bcftools index "${i}"
        done

        bcftools concat -a *.vcf.gz -O z -o !{output}

        get_nextflow_log.bash concat_vcf.log
        '''
    }
    //
}
