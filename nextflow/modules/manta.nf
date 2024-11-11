process MANTA {
    ext version: "1.6.0"
    conda { task.ext.conda }

    publishDir "$meta.out", mode: "copy"
    publishDir "$meta.log", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(tumor), path(normal)
    val(reference)
    val(is_exome)
    val(module_number)
    //

    output:
    tuple val(meta), path("*.vcf.gz"), emit: variants
    tuple val(meta.id), path("${module_number}-candidateSmallIndels.vcf.gz"), emit: indels
    path(out)
    path("*.log")
    //

    shell:
    out = "${module_number}-${meta.baseName}_MantaOut"
    check = file("${meta.out}/${out}")
    def exome_flag = is_exome ? " --exome " : ""
    if (check.exists()) {
        '''
        cp -r !{check}.name .
        cp !{meta.out}/*_Manta.vcf.gz .
        cp !{meta.log}/manta.log .
        '''
    } else {
        template 'manta.sh'
    }
    //
}
