include { FASTP } from "../modules/fastp.nf"
include { BWA } from "../modules/bwa.nf"
include { MARK_DUPLICATES } from "../modules/mark_duplicates.nf"
include { BQSR } from "../modules/bqsr.nf"
include { MUTECT2 } from "../modules/mutect2.nf"
include { MANTA } from "../modules/manta.nf"
include { STRELKA2 } from "../modules/strelka2.nf"
include { MSISENSORPRO } from "../modules/msisensorpro.nf"
include { CNVKIT } from "../modules/cnvkit.nf"
include { DELLY } from "../modules/delly.nf"
include { PICARD } from "../modules/picard.nf"

workflow "whole_exome" {

    main:
    def check_keys = params.ref.keySet() == ["genome",
                                             "homopolymers_microsatellites",
                                             "genome_copy_number",
                                             "genome_exclude",
                                             "known_variants",
                                             "panel_of_normals"].toSet()

    input = Channel.fromPath(params.input)
        .splitCsv(header: true)
        .map { it -> [["id": it.patient,
                "out": "${params.outdir}/${it.patient}_${it.source}",
                "type": it.source,
                "log": "${params.logdir}/${it.patient}_${it.source}",
                "RGLB": it.read_group_library, // Read group info used by GATK tools
                "RGPL": it.read_group_platform,
                "RGPU": it.read_group_platform_unit,
                "RGSM": it.read_group_sample_name
            ], [file(it.fastq_1), file(it.fastq_2)]] }

    /*
     * Preprocessing
     */

    FASTP(input, 1)
    BWA(FASTP.out.passed, params.ref.genome, 2)
    MARK_DUPLICATES(BWA.out.mapped, 3)
    BQSR(MARK_DUPLICATES.out.dedup, params.ref.genome, params.ref.known_variants, 4)

    // After preprocessing, branch data into tumor and normal then pair up by id
    output = BQSR.out.bam
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
    paired = tumors.join(normal)
    paired_no_id = paired.map { it -> [it[1], it[2], it[3]] }
    // paired_no_id is [meta, tumor, normal]

    /*
     * Variant calling
     */
    MUTECT2(paired, params.ref.genome, 4)
    MANTA(paired, params.ref.genome, true, 4)

    to_strelka = paired_no_id.join(MANTA.out.indels)
        .map { it -> [it[1], it[2], it[3], it[4]] }
    STRELKA(to_strelka, params.ref.genome, true, 4)

    // TODO need delly exclusion file
    DELLY(paired, params.ref.genome, params.ref.genome_exclude, 4)
    // TODO need cnvkit copy number file
    CNVKIT(paired, params.ref.genome_copy_number, 4)
    MSISENSORPRO(paired, params.ref.homopolymers_microsatellites, 4)


    /*
     * Variant annotation
     */



    /*
     * Metric collection
     */
}
