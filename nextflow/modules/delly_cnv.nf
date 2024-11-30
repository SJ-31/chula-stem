process DELLY_CNV {
    ext version: "1.3.1"
 
    publishDir "$meta.out", mode: "copy", saveAs: params.saveFn
    publishDir "$meta.log", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(normal), path(tumor), path(indices, arity: "2"), path(covfile), val(purity), val(ploidy)
    val(reference)
    val(mappability)
    val(module_number)
    //

    output:
    tuple val(with_caller), path(out)
    path(segmentation)
    path("*.log")
    //

    shell:
    out = params.getName(module_number, meta, "DellyCNV", ".vcf.gz")
    segmentation = params.getName(module_number, meta, "DellySegmentation", "vcf.gz")
    with_caller = meta + ["caller": "dellyCNV"]
    check = file("${meta.out}/${out}")
    check2 = file("${meta.out}/${segmentation}")
    args = task.ext.args.join(" ")
    if (check.exists() && check2.exists()) {
        '''
        ln -sr !{check} .
        ln -sr !{check2} .
        ln -sr !{meta.log}/dellyCNV.log .
        '''
    } else {
        template 'delly_cnv.bash'
    }
    //
}
