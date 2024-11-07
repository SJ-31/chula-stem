process MANTA {
    conda params.manta_env

    publishDir "$meta.out", mode: "copy"
    publishDir "$meta.log", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(tumor), path(normal)
    val(reference)
    val(is_exome)
    val(module_number)
    //

    output:
    tuple val(meta), path("*.vcf.gz"), emit: variants
    tuple val(meta.id), path("${module_number}-candidateSmallIndels.vcf.gz"), emit: indels
    path(out)
    path("*.log")
    //

    script:
    out = "${module_number}-${meta.baseName}_MantaOut"
    def check = file("${meta.out}/${out}")
    def exome_flag = is_exome ? " --exome " : ""
    if (check.exists()) {
        """
        cp $check.name .
        cp ${meta.log}/manta.log .
        """
    } else {
        """
        configManta.py \
            --normalBam $normal \
            --tumorBam $tumor \
            --referenceFasta $reference \
            $exome_flag \
            --runDir $out

        ./runWorkflow.py

        mv ${out}/variants/*.vcf.gz .
        for variant in *.vcf.gz; do
            mv \$variant "${module_number}-\${variant}"
        done

        bcftools annotate temp.vcf.gz TODO
        cp .command.out manta.log
        """
    }
    //
}
