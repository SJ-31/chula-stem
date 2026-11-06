process OPTITYPE {
    ext version: "1.5.0"

    input:
    // record(id: String), val(meta), path(reads)
    val data_type

    output:
    tuple val(meta), path(out_csv), emit: csv
    tuple val(meta), path(out_pdf), emit: pdf

    script:
    prefix = Utl.getName(meta, "optitype")
    out_csv = "${prefix}/"
    out_pdf = "${prefix}/"
    data_type_flag = "dna" == data_type ? "--dna" : "--rna"
    args = task.ext.args.join(" ")
    """
    optitype run -i "${reads[0]}" -i "${reads[1]}" \
        "${data_type_flag}" \
        -o . \
        --prefix "${prefix}" \
        ${args}
    """
}
