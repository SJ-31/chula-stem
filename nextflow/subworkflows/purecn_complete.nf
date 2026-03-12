include { MUTECT2 } from "../modules/mutect2.nf"
include { PURECN_BAIT_INTERVALS } from "../modules/purecn_bait_intervals.nf"
include { PURECN_COVERAGE as COVERAGE_NORMAL } from "../modules/purecn_coverage.nf"
include { PURECN_COVERAGE as COVERAGE_TUMOR } from "../modules/purecn_coverage.nf"
include { PURECN_NORMALDB } from "../modules/purecn_normaldb.nf"
include { PURECN_CALL } from "../modules/purecn_call.nf"


workflow PURECN_COMPLETE {
    //
    take:
    mutect2_unfiltered
    tumor_bam
    normal_bam
    panel_of_normals
    cohort_name
    module_number

    main:

    tumor_bam.dump(tag: "purecn_tumor")
    normal_bam.dump(tag: "purecn_normal")

    if (!params.ref.purecn_bait_intervals) {
        PURECN_BAIT_INTERVALS([["out": "${params.outdir}/PureCN_ref",
                                "log": "${params.outdir}/PureCN_ref"],
                               params.ref.genome],
                              params.ref.baits,
                              params.ref.mappability,
                              module_number)
        bait_intervals = PURECN_BAIT_INTERVALS.out.baits
    } else {
        bait_intervals = params.ref.purecn_bait_intervals
    }

    COVERAGE_TUMOR(tumor_bam, bait_intervals, false, module_number)

    if (!params.ref.purecn_normaldb) {
        COVERAGE_NORMAL(normal_bam, bait_intervals, false, module_number)
        normal_files = COVERAGE_NORMAL.out.loess.map({ it -> it[1] }).toList()
        def to_normaldb = channel.of(["filename": cohort_name,
                                      "out": "${params.outdir}/PureCN_ref",
                                      "log": "${params.outdir}/PureCN_ref"])
            .merge(normal_files) { meta, cov -> tuple(meta, tuple(cov)) }
        to_normaldb.dump(tag: "purecn_normaldb")
        PURECN_NORMALDB(to_normaldb, panel_of_normals, module_number)
        normaldb = PURECN_NORMALDB.out.db.first()
        mapping_bias = PURECN_NORMALDB.out.mapping_bias.first()
    } else {
        normaldb = params.ref.purecn_normaldb
        mapping_bias = params.ref.purecn_mapping_bias ?: ""
    }
    to_call = COVERAGE_TUMOR.out.loess.map({ it -> [it[0].id,
                                               ["type": "paired",
                                                "id": it[0].id,
                                                "filename": it[0].id],
                                               it[1]]})
        .join(Utl.getId(mutect2_unfiltered))
        .map({ it -> it[1..-1] })
        .map({ it ->
                [it[0] + ["out": "${params.outdir}/${it[0].id}/annotations",
                          "log": "${params.logdir}/${it[0].id}/annotations"]] +
                    it[1..-1] })

    snp_blacklist = params.ref.snp_blacklist ?: ""
    PURECN_CALL(to_call,
                normaldb,
                bait_intervals,
                mapping_bias,
                snp_blacklist,
                params.defaults.purity,
                params.defaults.ploidy,
                module_number + 1)

    emit:
    purity_ploidy = PURECN_CALL.out.purity_ploidy
}
