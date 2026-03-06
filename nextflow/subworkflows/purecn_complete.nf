include { MUTECT2 } from "../modules/mutect2.nf"
include { PURECN_BAIT_INTERVALS } from "../modules/purecn_bait_intervals.nf"
include { PURECN_COVERAGE } from "../modules/purecn_coverage.nf"
include { PURECN_NORMALDB } from "../modules/purecn_normaldb.nf"
include { PURECN_CALL } from "../modules/purecn_call.nf"


workflow PURECN_COMPLETE {
    //
    take:
    mutect2_unfiltered
    tumor_bam
    normal_bam
    module_number

    main:
    def cohort_name = params.cohort ? params.cohort : "cohort"

    if (params.ref.purecn_bait_intervals == null) {
        bait_intervals = PURECN_BAIT_INTERVALS(params.ref.genome,
                                               params.ref.baits,
                                               params.ref.mappability,
                                               module_number).out.baits.first()
    } else {
        bait_intervals = file(params.ref.purecn_bait_intervals) 
    }
    
    def tumor_cov = PURECN_COVERAGE(tumor_bam, bait_intervals, false, module_number)

    if (params.ref.purecn_normaldb == null) {
        def normal_cov = PURECN_COVERAGE(normal_bam, bait_intervals, false,
                                         module_number)
        def to_normaldb = channel.of(["filename": cohort_name,
                                      "out": "${params.outdir}/PureCN",
                                      "log": "${params.outdir}/PureCN"])
            .merge(normal_cov.out.loess.toList()) { meta, cov ->
                tuple(meta, tuple(cov)) }
        normaldb = PURECN_NORMALDB(to_normaldb, params.ref.panel_of_normals,
                                   module_number).out.db.first() 
    } else {
        normaldb = file(params.ref.purecn_normaldb)       
    }
    to_call = tumor_cov.out.loess.map({ it -> [it[0].id,
                                               ["type": "paired",
                                                "id": it[0].id,
                                                "filename": it[0].id],
                                               it[1]]})
        .join(Utl.getId(mutect2_unfiltered.out.variants))
        .map({ it -> it[1..-1] })
    
    PURECN_CALL(to_call, normaldb, module_number + 1)
}
