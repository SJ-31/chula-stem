include { PURECN_BAIT_INTERVALS } from "../modules/purecn_bait_intervals.nf"
include { PURECN_COVERAGE as COVERAGE_NORMAL } from "../modules/purecn_coverage.nf"
include { PURECN_COVERAGE as COVERAGE_TUMOR } from "../modules/purecn_coverage.nf"
include { PURECN_NORMALDB } from "../modules/purecn_normaldb.nf"
include { PURECN_CALL } from "../modules/purecn_call.nf"

nextflow.enable.types = true

record Sample {
    id: String
    mutect2_vcf: Path
    tumor_bam: Path
    normal_bam: Path
    tumor_bam_index: Path
    normal_bam_index: Path
}

workflow PURECN_COMPLETE {
    take:
    samples: Channel<Sample>
    panel_of_normals: String

    main:
    if (!params.ref.purecn_bait_intervals) {
        bait_intervals = PURECN_BAIT_INTERVALS(
            params.ref.genome,
            params.ref.baits,
            params.ref.mappability,
        ).toList().map { v -> v.first() }
        // NOTE: [2026-11-07 Sat] as of Nextflow version 26.10.5
        // there doesn't seem to be a better idiom for generating a value
        // channel. The type checker will complain about the above, but
        // it runs fine
    }
    else {
        bait_intervals = channel.value(params.ref.purecn_bait_intervals)
    }
    tumor_bam = samples.map { v -> record(id: v.id, bam_or_cov: v.tumor_bam, index: v.tumor_bam_index) }
    normal_bam = samples.map { v -> record(id: v.id, bam_or_cov: v.normal_bam, index: v.normal_bam_index) }


    if (!params.ref.purecn_normaldb && !params.ref.purecn_mapping_bias) {
        to_normaldb = COVERAGE_NORMAL(normal_bam, bait_intervals, false).map { val -> val.loess }.collect()

        db_ch = PURECN_NORMALDB(to_normaldb, panel_of_normals, params.genome_build).toList().map { b -> b.first() }
        normaldb = db_ch.map { v -> v.db }
        mapping_bias = db_ch.map { v -> v.mapping_bias }
    }
    else {
        normaldb = params.ref.purecn_normaldb
        mapping_bias = params.ref.purecn_mapping_bias
    }
    coverage_tumor = COVERAGE_TUMOR(tumor_bam, bait_intervals, false)

    to_call = coverage_tumor.map { v -> record(id: v.id, coverage: v.coverage) }.join(samples.map { v -> record(id: v.id, mutect2_vcf: v.mutect2_vcf) }, by: "id")

    purecn_call = PURECN_CALL(
        to_call,
        normaldb,
        bait_intervals,
        mapping_bias,
        params.ref.snp_blacklist,
        params.defaults.purity,
        params.defaults.ploidy,
    )

    emit:
    bait_intervals = bait_intervals
    call = purecn_call
    coverage = coverage_tumor
}
