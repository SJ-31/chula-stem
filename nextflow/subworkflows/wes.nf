workflow "whole_exome" {

    main:
    Channel.fromPath(params.input)
        .splitCsv(header: true)
        .map { it -> [it.patient, [it.fastq_1, it.fastq_2]] }
        .set { input }

    println params.name
    input.view()

}
