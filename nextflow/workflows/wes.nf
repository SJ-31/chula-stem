include { FASTP } from "../modules/fastp.nf"
include { BWA } from "../modules/bwa.nf"
include { MARK_DUPLICATES } from "../modules/mark_duplicates.nf"
include { BQSR } from "../modules/bqsr.nf"
include { MUTECT2 } from "../modules/mutect2.nf"
include { MANTA } from "../modules/manta.nf"
include { STRELKA2 } from "../modules/strelka2.nf"
include { MSISENSORPRO } from "../modules/msisensorpro.nf"
include { CNVKIT } from "../modules/cnvkit.nf"
include { DELLY_SV } from "../modules/delly_sv.nf"
include { PICARD } from "../modules/picard.nf"
include { SNPEFF } from "../modules/snpeff.nf"
include { SNPSIFT } from "../modules/snpsift.nf"
include { MERGE_VCFS as MERGE_A } from "../modules/merge_vcfs.nf"
include { MERGE_VCFS as MERGE_B } from "../modules/merge_vcfs.nf"
include { DELLY_COV } from "../modules/delly_cov.nf"
include { VEP } from "../modules/vep.nf"
include { SAMTOOLS_INDEX } from "../modules/samtools_index.nf"

workflow "whole_exome" {

    main:
    input = Channel.fromPath(params.input)
        .splitCsv(header: true)
        .map { [["id": it.patient,
                "out": "${params.outdir}/${it.patient}_${it.source}",
                "type": it.source,
                "log": "${params.logdir}/${it.patient}_${it.source}",
                "filename": "${it.patient}_${it.source}", // Output filename
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
    MARK_DUPLICATES(BWA.out.mapped, params.ref.genome, 3)
    BQSR(MARK_DUPLICATES.out.dedup, params.ref.genome, params.ref.known_variants, 4)

    // After preprocessing, branch data into tumor and normal then pair up by id and join with
    //  index
    output = BQSR.out.bam
    SAMTOOLS_INDEX(output) // Required by certain callers
    branched = output.branch { meta, files ->
        tumor: meta.type == "cancer" || meta.type == "tumor"
        normal: meta.type == "normal"
    }
    normals = branched.normal.map { [it[0].id,
                                         ["type": "paired",
                                          "id": it[0].id,
                                          "out": "${params.outdir}/${it[0].id}_paired",
                                          "RGSM_normal": it[0].RGSM,
                                          "filename": it[0].id,
                                          "log": "${params.logdir}/${it[0].id}_paired"],
                                         it[1]] }
    tumors = branched.tumor.map { [it[0].id, it[1]] }
    indices = SAMTOOLS_INDEX.out.index.map({ [it[0].id, it[1]] }).groupTuple()
    paired = normals.join(tumors).join(indices) // Order is important for the first two
    paired_no_id = paired.map { it[1..-1] }
    // paired_no_id is [meta, normal, tumor, indices]

    // DELLY_COV(branched.tumor, params.ref.genome, params.ref.mappability, 3)
    /*
     * Variant calling
     */
    def getId = { [it[0].id] + it[1] }

    MUTECT2(paired_no_id, params.ref.genome, 4)
    MANTA(paired_no_id, params.ref.genome, 4)

    to_strelka = paired.join(MANTA.out.indels)
        .map { it[1..-1] }
    STRELKA2(to_strelka, params.ref.genome, 4)

    // to_delly_cnv = paired.join(DELLY_COV.out.cov).map {it[1..-1]}
    // DELLY_CNV(to_delly_cnv, params.ref.genome, params.ref.mappability, 4)

    DELLY_SV(paired_no_id, params.ref.genome, params.ref.delly_exclude, 4)
    // TODO need cnvkit copy number file
    // CNVKIT(paired_no_id, params.ref.cnvkit_copy_number, 4)
    // MSISENSORPRO(paired_no_id, params.ref.homopolymers_microsatellites, 4)
    //
    // all_variants = MUTECT2.out.variants.join(
    //     MANTA.out.somatic.map(getId),
    //     STRELKA2.out.somatic.map(getId),
    //     DELLY_SV.out.variants.map(getId),
    //     DELLY_CNV.out.variants.map(getId),
    //     // TODO: This is incomplete cause you don't know the output of the other callers
    // ).map { it[1..-1] } MERGE_A(all_variants, 5)

    /*
     * Variant annotation
     */
    // SNPEFF(CONCAT_VCF_1.out.vcf, params.ref.snpEff_db, true, 6)
    // VEP(all_variants.out.vcf, params.ref.genome, 6)
    // TODO: need to transfer annotations between them
    //      Check if merge1 works for this
    // annotated = SNPEFF.out.vcf.join(VEP.map(getId))
    // MERGE_B(annotated, 6)

    /*
     * Metric collection
     */


}
