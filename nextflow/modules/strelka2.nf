process STRELKA2 {
    ext version: "2.9.10"

    publishDir "$meta.out", mode: "copy"
    publishDir "$meta.log", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(tumor), path(normal), path(manta_indels)
    val(reference)
    val(module_number)
    //

    output:
    tuple val(meta), path("*.vcf.gz"), emit: variants
    path(out)
    path("*.log")
    //

    script:
    out = "${module_number}-${meta.id}_StrelkaOut"
    check = file("${meta.out}/${out}")
    def args = task.ext.args.join(" ")
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
