process STRELKA2 {
    ext version: "2.9.10"

    publishDir "$meta.out", mode: "copy", saveAs: { x -> x ==~ /.*\.log/ ? null : x }
    publishDir "$meta.log", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(normal), path(tumor), path(indices, arity: "2"), path(manta_indels)
    val(reference)
    val(module_number)
    //

    output:
    tuple val(meta), path("*.vcf.gz"), emit: variants
    path(out)
    path("*.log")
    //

    shell:
    out = "${module_number}-${meta.id}_StrelkaOut"
    check = file("${meta.out}/${out}")
    args = task.ext.args.join(" ")
    if (check.exists()) {
        '''
        ln -sr !{check} .
        ln -sr !{meta.out}/*_Strelka.vcf.gz .
        ln -sr !{meta.log}/strelka.log .
        '''
    } else {
        template 'strelka.sh'
    }
    //
}
