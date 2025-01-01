process CONCAT_VCF {
    ext version: "1.21"
    // Concatenate groups vcf files of the same samples

    publishDir "$meta.out", mode: "copy", saveAs: params.saveFn
    publishDir "$meta.log", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(vcfs, arity: "2..*")
    val(reference)
    val(module_number)
    //

    output:
    tuple val(meta), path(output), emit: vcf
    path("*.log")
    //

    script:
    output = Utl.getName(module_number, meta, "All", "vcf.gz")
    check = file("${meta.out}/${output}")
    if (check.exists()) {
        """
        ln -sr ${check} .
        ln -sr ${meta.log}/concat_vcf.log .
        """
    } else {
        """
        i=0
        for v in *.vcf.gz; do
            name="\${i}.bcf.gz"
            bcftools annotate -x 'INFO/HOMLEN,INFO/SVLEN,FORMAT/SR' "\${v}" -O b > "\${name}"
            bcftools index "\${name}"
            i=\$((i+1))
        done

        bcftools concat -a *.bcf.gz -O z -o ${output}

        get_nextflow_log.bash concat_vcf.log
        """
    }
    //
}
