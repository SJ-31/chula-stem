nextflow.enable.types = true

process CONCAT_VCF {
    input:
    // Concatenate groups vcf files of the same samples
    record(
        id: String,
        vcfs: Bag<Path>
    )
    prefix: String

    output:
    record(
        id: id,
        concat: output,
    )

    topic:
    record(id: id, process: "concat_vcf", file: file("concat_vcf.log")) >> "logs"

    script:
    //
    output = "${id}-${prefix}concat.vcf.gz"
    """
    i=0
    for v in *.vcf.gz; do
        name="\${i}.bcf.gz"
        if [[ -s "\${v}" ]]; then
            bcftools annotate -x 'INFO/HOMLEN,INFO/SVLEN,FORMAT/SR' "\${v}" -O b > "\${name}"
            bcftools index "\${name}"
            i=\$((i+1))
         fi
    done

    bcftools concat -a *.bcf.gz -O z -o ${output}

    get_nextflow_log.bash concat_vcf.log
    """
}
