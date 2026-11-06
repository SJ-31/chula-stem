nextflow.enable.types = true

process OPTITYPE {
    
    input:
    record(id: String,
           fastq_1: Path,
           fastq_2: Path
    )
    data_type: String

    output:
    record(id: id,
           optitype_csv: out_csv,
           optitype_pdf: out_pdf
    )

    script:
    prefix = "${id}-optitype"
    out_csv = "${prefix}/"
    out_pdf = "${prefix}/"
    data_type_flag = "dna" == data_type ? "--dna" : "--rna"
    args = task.ext.args.join(" ")
    """
    optitype run -i "${fastq_1}" -i "${fastq_2}" \
        "${data_type_flag}" \
        -o . \
        --prefix "${prefix}" \
        ${args}
    """
}
