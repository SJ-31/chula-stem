process DELLY_SV {
    ext version: "1.3.1"
 
    publishDir "$meta.out", mode: "copy", saveAs: { x -> x ==~ /.*\.log/ ? null : x }
    publishDir "$meta.log", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(normal), path(tumor), path(indices, arity: "2")
    val(reference)
    val(exclude)
    val(module_number)
    //

    output:
    tuple val(meta), path(out), emit: variants
    path("*.log")
    //

    shell:
    out = "${module_number}-${meta.filename}-DellySV.vcf.gz"
    // the output of this is yet
    check = file("${meta.out}/${out}")
    args = task.ext.args.join(" ")
    exclude_flag = exclude == "" ? "" : "-x ${exclude}"
    if (check.exists()) {
        '''
        ln -sr !{check} .
        ln -sr !{meta.log}/dellySV.log .
        '''
    } else {
        template 'delly.bash'
    }
    //
}
