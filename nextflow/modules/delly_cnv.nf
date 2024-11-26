process DELLY_CNV {
    ext version: "1.3.1"
 
    publishDir "$meta.out", mode: "copy", saveAs: params.saveFn
    publishDir "$meta.log", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(normal), path(tumor), path(indices, arity: "2"), path(covfile)
    val(reference)
    val(mappability)
    val(module_number)
    //

    output:
    tuple val(meta), path(out)
    path("*.log")
    //

    shell:
    out = "${module_number}-${meta.filename}-DellyCNV.vcf.gz"
    // the output of this is yet
    check = file("${meta.out}/${out}")
    args = task.ext.args.join(" ")
    if (check.exists()) {
        '''
        ln -sr !{check} .
        ln -sr !{meta.log}/dellyCNV.log .
        '''
    } else {
        template 'delly_cnv.bash'
    }
    //
}
