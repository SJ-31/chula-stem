process DNDSCV {
    ext version: "0.1.0"

    publishDir "$meta.out", mode: "copy"
    publishDir "$meta.log", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(vcf)
    // The input format is a tsv file with five columns: sampleID, chromosome, position,
    //      reference base and mutant base
    // We generate this directly from the vcf (dndscv supports snps and indels)

    output:
    path("*.log")
    //

    script:
    output =
        check = file("${meta.out}/output")
    def args = task.ext.args.join(" ")
    if (check.exists()) {
        """

        cp ${meta.log}/ .
        """
    } else {
        """
        echo -e "chr\tpos\tref\tmut" > table.tsv
        bcftools query \\
            -i "SNP=1 || INS=1 || DEL=1" $vcf \\
            -f "%CHROM\t%POS0\t%REF\t%ALT" \\
            >> table.tsv

        cp .command.out .log
        """
    }
    //
}
