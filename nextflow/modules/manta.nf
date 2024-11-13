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
    tuple val(meta), path("${module_number}-somaticSV_Manta.vcf.gz"), emit: somatic
    tuple val(meta.id), path("${module_number}-candidateSmallIndels_Manta.vcf.gz"), emit: indels
    path(out)
    path("*.log")
    // Produces 4 files
    // - diploidSV: variants called under a diploid model, i.e. the normal sample only
    // - somaticSV: variants called under somatic model, considering paired vs. normal
    // - candidateSV: unscored SV and indel candidates
    // - candidateSmallIndels: subset of candidateSV with only small indels (less than 50 by default)

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
