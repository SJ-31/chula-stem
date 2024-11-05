include { FASTP } from "../modules/fastp.nf"

workflow "whole_exome" {

    main:
    Channel.fromPath(params.input)
        .splitCsv(header: true)
        .map { it -> [["sample": it.patient,
                "out": "${params.outdir}/${it.patient}",
                "log": "${params.logdir}/${it.patient}",
            ], [file(it.fastq_1), file(it.fastq_2)]] }
        .set { input }

   FASTP(input, 1)
}
