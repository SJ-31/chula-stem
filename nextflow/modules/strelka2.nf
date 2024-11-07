process STRELKA {
    // Version 2.9.10
    conda params.strelka_env

    publishDir "$meta.out", mode: "copy"
    publishDir "$meta.log", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(tumor), path(normal), path(manta_indels)
    val(reference)
    val(is_exome)
    val(module_number)
    //

    output:
    tuple val(meta), path("*.vcf.gz"), emit: variants
    path(out)
    path("*.log")
    //

    script:
    out = "${module_number}-${meta.baseName}_StrelkaOut"
    def check = file("${meta.out}/${out}")
    def exome_flag = is_exome ? " --exome " : ""
    if (check.exists()) {
        '''
        cp -r !{check}.name .
        cp !{meta.out}/*_Strelka.vcf.gz .
        cp !{meta.log}/strelka.log .
        '''
    } else {
        template 'strelka.sh'
    }
    //
}
