process SNPEFF {
    ext version: "5.2e"
    
    publishDir "$meta.out", mode: "copy", saveAs: params.saveFn
    publishDir "$meta.log", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(vcf)
    val(reference)
    val(is_somatic)
    val(module_number)
    //

    output:
    tuple val(meta), path("${output}.gz"), emit: vcf
    tuple val(meta), path("*.html"), path("*.txt"), emit: report
    path("*.log")
    //

    shell:
    output = params.getName(module_number, meta, "snpeFF", "vcf")
    report = params.getName(module_number, meta, "snpEff_summary", "html")
    genes_file = params.getName(module_number, meta, "snpEff_genes", "txt")
    check = file("$meta.out/${output}.gz")
    args = task.ext.args.join(" ")
    if (check.exists()) {
        '''
        ln -sr !{check} .
        ln -sr "!{meta.out}/!{report}" .
        ln -sr "!{meta.out}/!{genes_file}" .
        ln -sr "!{meta.log}/snpEff.log" .
        '''
    } else {
        template 'snpeff.sh'
    }
    //
}
