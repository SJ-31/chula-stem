process MANTA {
    ext version: "1.6.0"

    publishDir "$meta.out", mode: "copy"
    publishDir "$meta.log", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(normal), path(tumor), path(indices, arity: "2")
    val(reference)
    val(module_number)
    //

    output:
    tuple val(meta), path("*.vcf.gz"), emit: variants
    tuple val(meta.id), path("${module_number}-candidateSmallIndels.vcf.gz"), emit: indels
    path(out)
    path("*.log")
    //

    shell:
    out = "${module_number}-${meta.id}_MantaOut"
    check = file("${meta.out}/${out}")
    args = task.ext.args.join(" ")
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
