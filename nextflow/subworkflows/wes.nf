import groovy.toml.TomlSlurper

include { FASTP } from "../modules/fastp.nf"
include { BWA } from "../modules/bwa.nf"

workflow "whole_exome" {

    main:
    references = new TomlSlurper().parse(file(params.references))
    input = Channel.fromPath(params.input)
        .splitCsv(header: true)
        .map { it -> [["id": it.patient,
                "out": "${params.outdir}/${it.patient}_${it.source}",
                "type": it.source,
                "log": "${params.logdir}/${it.patient}_${it.source}",
            ], [file(it.fastq_1), file(it.fastq_2)]] }

   FASTP(input, 1)
   BWA(FASTP.out.passed, references.genome, 2)


    // TODO: after preprocessing, branch data into tumor and normal
    // then to pair up tumour-normal samples
    s1 = [["id": 1, "type": "tumor", "out": "foo", "log": "foo"], "bar"]
    s2 = [["id": 1, "type": "normal", "out": "foo", "log": "foo"], "bar"]
    output = Channel.from([s1, s2])
    branched = output.branch { it ->
        tumor: it[0].type == "tumor"
        normal: it[0].type == "normal"
    }
    tumors = branched.tumor.map { it -> [it[0].id,
                                         ["type": "paired", "id": it[0].id,
                                         "out": "${params.outdir}/${it[0].id}_paired",
                                         "log": "${params.logdir}/${it[0].id}_paired"],
                                         it[1]] }
    normal = branched.normal.map { it -> [it[0].id, it[1]] }
    paired = tumors.join(normal).map { it -> [it[1], it[2], it[3]] }

}
