process SNPEFF {
    ext version: "5.2e"
    
    publishDir "$meta.out", mode: "copy", saveAs: { x -> x ==~ /.*\.log/ ? null : x }
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
    output = "${module_number}-${meta.id}_snpEff.vcf"
    report = "${module_number}-${meta.id}_snpEff_summary.html"
    genes_file = "${module_number}-${meta.id}_snpEff_genes.txt"
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
