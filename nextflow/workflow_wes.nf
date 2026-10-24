include { MUTECT2_COMPLETE } from "./subworkflows/mutect2_complete.nf"
include { PREPROCESS_FASTQ } from "./subworkflows/preprocess_fastq.nf"
include { PURECN_COMPLETE } from "./subworkflows/purecn_complete.nf"
include { MANTA } from "./modules/manta.nf"
include { STRELKA2 } from "./modules/strelka2.nf"
include { MSISENSORPRO } from "./modules/msisensorpro.nf"
include { MSISENSORPRO_COLLECT } from "./modules/msisensorpro.nf"
include { CNVKIT } from "./modules/cnvkit.nf"
include { CNVKIT_PREP } from "./modules/cnvkit_prep.nf"
include { MOSDEPTH } from "./modules/mosdepth.nf"
include { REPORT } from "./modules/report.nf"
include { BCFTOOLS_STATS } from "./modules/bcftools_stats.nf"
include { MULTIQC } from "./modules/multiqc.nf"
include { DELLY_SV } from "./modules/delly_sv.nf"
include { PICARD } from "./modules/picard.nf"
include { CLASSIFY_CNV_FORMAT } from "./modules/classify_cnv_format.nf"
include { OCTOPUS } from "./modules/octopus.nf"
include { CLAIRS } from "./modules/clairs.nf"
include { CLASSIFY_CNV } from "./modules/classify_cnv.nf"
include { CONCAT_VCF as CONCAT_SMALL_1 } from "./modules/concat_vcf.nf"
include { CONCAT_VCF as CONCAT_SMALL_2 } from "./modules/concat_vcf.nf"
include { CONCAT_VCF as CONCAT_SV } from "./modules/concat_vcf.nf"
include { CALLSET_QC as QC_SMALL } from "./modules/callset_qc.nf"
include { CALLSET_QC_TSV } from "./modules/callset_qc_tsv.nf"
include { CROSS_REFERENCE as CROSS_REFERENCE_MSI } from "./modules/cross_reference.nf"
include { CROSS_REFERENCE as CROSS_REFERENCE_CNV } from "./modules/cross_reference.nf"
include { STANDARDIZE_VCF } from "./modules/standardize_vcf.nf"
include { MUSE2 } from "./modules/muse2.nf"
include { GRIDSS } from "./modules/gridss.nf"
include { SIGPROFILERASSIGNMENT } from "./modules/sigprofilerassignment.nf"
include { GET_THERAPY_CACHE } from "./modules/get_therapy_cache.nf"
include { SIGPROFILERASSIGNMENT_COLLECT } from "./modules/sigprofilerassignment.nf"
include { VEP } from "./modules/vep.nf"
include { PANEL_OF_NORMALS_FROM_BAM } from "./subworkflows/panel_of_normals_from_bam.nf"


def toConcat(suffix, val) {
    [val[0] + ["suffix": suffix]] + [val[1..-1]]
}

params {
    tumor_only = false
}

workflow {

    main:
    def cohort_name = params.cohort ?: "cohort"
    def branchSources = branchCriteria { it ->
        tumor: it[0].type == "cancer" || it[0].type == "tumor"
        normal: it[0].type == "normal"
    }
    def cohortTopLevel = { it -> [["filename": cohort_name], it] }
    /*
     * Preprocessing
     */

    PREPROCESS_FASTQ(params.input, "wes")
    // After preprocessing, branch data into tumor and normal then pair up by id and join with  index

    branched = PREPROCESS_FASTQ.out.bam.branch(branchSources)
    
    normals = Utl.getId(
        branched.normal.map { val ->
            [["type": "paired", "id": val[0].id, "RGSM_normal": val[0].RGSM, "filename": val[0].id], val[1]]
        }
    )
    tumors = branched.tumor.map { val -> [["RGSM_tumor": val[0].RGSM], val[1]] }
    
    indices = Utl.getId(PREPROCESS_FASTQ.out.bam_index).groupTuple()
    paired = normals
        .join(tumors)
        .map { it -> [it[0]] + [it[1] + it[3]] + [it[2]] + [it[4]] } // Merge the maps
        .join(indices)
    // Order is important for the first two
    paired_no_id = Utl.delId(paired)
    // paired_no_id is [meta, normal, tumor, [normal_index, tumor_index]]
    // the order of indices doesn't matter, they just need to be present

    /*
     * (Optional) Panel of normals
     */
    if (params.create_pon) {
        to_pon_bam = branched.normal.map { it ->
            [it[0] + ["out": "${params.outdir}/for_panel_of_normals"]] + it[1..-1]
        }
        to_pon_bam_indices = PREPROCESS_FASTQ.out.bam_index.filter { v -> v[0].type == "normal" }

        to_pon_bam.dump(tag: "pon_bam")
        to_pon_bam_indices.dump(tag: "pon_bam_indices")
        PANEL_OF_NORMALS_FROM_BAM(
            to_pon_bam,
            to_pon_bam_indices,
            cohort_name,
            channel.empty(),
            false,
        )
        panel_of_normals = PANEL_OF_NORMALS_FROM_BAM.out.pon_ws.map { it -> it[1] }
    }
    else {
        panel_of_normals = params.ref.panel_of_normals ?: ""
    }

    /*
     * Variant calling
     */

    // Structural variants

    MANTA(paired_no_id, params.ref.genome, params.ref.targets)
    DELLY_SV(paired_no_id, params.ref.genome, params.ref.delly_exclude)

    MSISENSORPRO(
        paired_no_id,
        params.ref.homopolymers_microsatellites,
        "exome",
        params.ref.genome_gff,
    )

    MSISENSORPRO_COLLECT(MSISENSORPRO.out.summary.collect().map(cohortTopLevel))
    // TODO: think about trying out https://github.com/walaj/svaba
    // GRIDSS(paired_no_id, params.ref.genome, params.ref.genome_blacklist, 5)

    // Small variants
    to_mutect = paired_no_id.map { it ->
        [it[0] + ["out": "${it[0].out}/5-Mutect2"]] + it[1..-1]
    }
    MUTECT2_COMPLETE(to_mutect, panel_of_normals)
    to_strelka = Utl.delId(Utl.joinFirst(paired, [MANTA.out.indels]))
    STRELKA2(to_strelka, params.ref.genome, params.ref.targets)
    MUSE2(paired_no_id, params.ref.genome, params.ref.dbsnp, "exome")

    // Combine variants by type
    small_variants_to_geno = Utl
        .joinFirst(
            MUTECT2_COMPLETE.out.filtered,
            [STRELKA2.out.variants.map { it -> [it[0]] + it[1] }, MUSE2.out.variants],
        )
        .map { it -> toConcat("Concat_to_oct", it) }
    CONCAT_SMALL_1(small_variants_to_geno)

    // Octopus and Clairs uses previous variants to aid calling
    to_geno = Utl.joinFirst(paired_no_id, [CONCAT_SMALL_1.out.vcf])
    OCTOPUS(to_geno, params.ref.genome, params.ref.targets, false)
    CLAIRS(to_geno, params.ref.genome, params.ref.targets)

    to_concat_small_2 = Utl
        .joinFirst(
            CONCAT_SMALL_1.out.vcf,
            [OCTOPUS.out.variants, CLAIRS.out.variants],
        )
        .map { it -> [it[0]] + [it[1..-1]] }
    CONCAT_SMALL_2(to_concat_small_2)
    structural_variants = Utl
        .joinFirst(
            MANTA.out.somatic,
            [
                DELLY_SV.out.variants
                // GRIDSS.out.variants
            ],
        )
        .map { it -> toConcat("SV_all", it) }
    CONCAT_SV(structural_variants)
    to_std_small = Utl.joinFirst(
        CONCAT_SMALL_2.out.vcf,
        [branched.normal, branched.tumor],
    )

    /*
     * QC, Re-compute AD, DP and VAF for small variants (some callers do not compute this innately)
     */

    STANDARDIZE_VCF(Utl.addSuffix(to_std_small, "Small_std"), params.ref.genome)

    small_all = Utl.delSuffix(STANDARDIZE_VCF.out.vcf)
    sv_all = Utl.delSuffix(CONCAT_SV.out.vcf)

    QC_SMALL(
        Utl.addSuffix(small_all, "Small_high_conf"),
        params.small_qc,
        panel_of_normals,
    )

    /*
     * Copy number aberration
     */

    if (!params.ref.cnvkit_reference) {
        collected_normals = normals
            .map { it -> it[2] }
            .mix(PREPROCESS_FASTQ.out.bam_index.map { it -> it[1] })
            .toList()
        CNVKIT_PREP(
            channel.of(
                ["filename": cohort_name, "out": "${params.outdir}/cnvkit_reference", "log": "${params.outdir}/cnvkit_reference"]
            ).merge(collected_normals) { meta, bams -> tuple(meta, tuple(bams)) },
            params.ref.genome,
            params.ref.baits_unzipped,
            params.ref.genome_blacklist,
            true,
            "hybrid",
        )
        cnvkit_reference = CNVKIT_PREP.out.reference.first()
        cnvkit_autobin = CNVKIT_PREP.out.autobin.first()
    }
    else {
        cnvkit_reference = params.ref.cnvkit_reference
        cnvkit_autobin = params.ref.cnvkit_autobin
    }

    branched_with_index = PREPROCESS_FASTQ.out.bam_with_index.branch(branchSources)
    fastp = PREPROCESS_FASTQ.out.fastp_json.branch(branchSources)
    trimmed = PREPROCESS_FASTQ.out.trimmed.branch(branchSources)

    PURECN_COMPLETE(
        MUTECT2_COMPLETE.out.unfiltered,
        branched_with_index.tumor,
        branched_with_index.normal,
        panel_of_normals,
        cohort_name,
    )

    to_cnvkit = Utl.delId(

        paired.map { it -> it[0..1] + [it[3]] }.join(Utl.getId(MUTECT2_COMPLETE.out.filtered)).join(Utl.getId(PURECN_COMPLETE.out.purity_ploidy))
    )

    CNVKIT(to_cnvkit, cnvkit_reference, "hybrid")

    CLASSIFY_CNV_FORMAT(CNVKIT.out.cns)
    cnv_bed = CLASSIFY_CNV_FORMAT.out.bed
        .collectFile(
            keepHeader: true,
            skip: 1,
        ) { meta, file -> ["5-${meta.id}-ClassifyCNV_all.bed", file] }
        .map { it ->
            def id = (it.baseName =~ /[0-9]+-(.*)-ClassifyCNV.*/)[0][1]
            [["filename": id, "id": id], it]
        }

    /*
     * Variant annotation
     */
    to_vep_small = Utl.modifyMeta(
        small_all,
        ["suffix": "VEP_small", "variant_class": "small", "qc": params.small_qc],
    )
    to_vep_sv = Utl.modifyMeta(
        sv_all,
        ["suffix": "VEP_SV", "variant_class": "sv", "qc": params.sv_qc],
    )

    VEP(to_vep_small.mix(to_vep_sv), params.ref.genome, panel_of_normals)

    to_qc_tsv = Utl.joinFirst(
        VEP.out.tsv,
        [MSISENSORPRO.out.tsv],
    )

    CALLSET_QC_TSV(to_qc_tsv, "")

    SIGPROFILERASSIGNMENT(
        Utl.delSuffix(QC_SMALL.out.vcf),
        true,
        "${params.configdir}/excluded_signatures.txt",
    )
    SIGPROFILERASSIGNMENT_COLLECT(
        SIGPROFILERASSIGNMENT.out.activities.map { it -> it[1] }.collect().map(cohortTopLevel)
    )
    CLASSIFY_CNV(cnv_bed)

    // Cross reference regions
    CROSS_REFERENCE_CNV(
        Utl.addSuffix(CLASSIFY_CNV.out.tsv, "CR_CNV"),
        "CNV",
        params.ref.clingen_dosage,
        params.ref.cnv_reference,
    )

    to_cr_msi = Utl.modifyMeta(MSISENSORPRO.out.tsv, [suffix: "CR_MSI"])
    CROSS_REFERENCE_MSI(
        to_cr_msi,
        "MSI",
        params.ref.clingen_gene,
        params.ref.msi_reference,
    )

    /*
     * Metric collection
     */

    to_metrics = PREPROCESS_FASTQ.out.bam_with_index

    PICARD(
        to_metrics,
        "hs",
        params.ref.genome,
        params.ref.targets_il,
        params.ref.baits_il,
        "",
        "",
    )
    MOSDEPTH(to_metrics, params.ref.targets)

    to_bcftools_stats = Utl
        .addSuffix(small_all, "Bcftools_stats_small")
        .mix(Utl.addSuffix(sv_all, "Bcftools_stats_SV"))

    BCFTOOLS_STATS(to_bcftools_stats, params.ref.targets)

    to_multiqc = PREPROCESS_FASTQ.out.fastp_json
        .mix(
            VEP.out.report,
            MOSDEPTH.out.dist,
            PICARD.out.metrics,
            BCFTOOLS_STATS.out.py,
        )
        .flatten()
        .collect()
        .map(cohortTopLevel)
    MULTIQC(to_multiqc)

    // Combine channels for report
    vep_grouped = Utl
        .getId(CALLSET_QC_TSV.out.tsv)
        .groupTuple(sort: { it -> (it.baseName =~ /VEP_small/) ? 1 : -1 })
    vep_to_report = vep_grouped.map { it -> [it[0]] + [vep_small: it[1][1], vep_sv: it[1][0]] }

    others = Utl
        .joinFirst(
            SIGPROFILERASSIGNMENT.out.activities,
            [CNVKIT.out.cns, CNVKIT.out.cnr, CROSS_REFERENCE_CNV.out.tsv, CROSS_REFERENCE_MSI.out.tsv],
        )
        .map { it ->
            [it[0], [sigprofiler: it[1], cnvkit_cns: it[2], cnvkit_cnr: it[3], classify_cnv: it[4], msisensor_pro: it[5], cosmic_reference: params.ref.cosmic_signatures]]
        }

    to_report = Utl.delId(Utl.getId(others, true).join(vep_to_report)).map { it -> [it[0], it[1] + it[2]] }

    if (params.cache_lookup) {
        GET_THERAPY_CACHE(
            vep_to_report.map { val -> val[1].vep_small }.collect(),
            vep_to_report.map { val -> val[1].vep_sv }.collect(),
            params.ref.therapy_data_cache,
        )
        therapy_data_cache = GET_THERAPY_CACHE.out.cache
    }
    else {
        therapy_data_cache = channel.of(params.ref.therapy_data_cache)
    }

    REPORT(to_report, "variant_calling")

    to_publish = Utl
        .joinFirst(
            REPORT.out.report,
            [
                branched.normal,
                branched.tumor,
                trimmed.normal,
                trimmed.tumor,
                fastp.normal,
                fastp.tumor,
                // Small variant callers
                MUTECT2_COMPLETE.out.unfiltered,
                MUTECT2_COMPLETE.out.filtered,
                MUTECT2_COMPLETE.out.stats,
                MUTECT2_COMPLETE.out.contamination,
                MUTECT2_COMPLETE.out.read_orientation,
                MUTECT2_COMPLETE.out.pileup,
                MUTECT2_COMPLETE.out.segmentation,
                STRELKA2.out.variants,
                STRELKA2.out.dir,
                MUSE2.out.variants,
                MUSE2.out.txt,
                OCTOPUS.out.variants,
                OCTOPUS.out.all,
                CLAIRS.out.variants,
                // Other callers
                MANTA.out.variants,
                MANTA.out.indels,
                MANTA.out.somatic,
                MANTA.out.dir,
                DELLY_SV.out.variants,
                MSISENSORPRO.out.all,
                PURECN_COMPLETE.out.call,
                PURECN_COMPLETE.out.loess,
                PURECN_COMPLETE.out.cov,
                CNVKIT.out.all,
                SIGPROFILERASSIGNMENT.out.all,
                // Annotated calls
                CROSS_REFERENCE_CNV.out.tsv,
                CROSS_REFERENCE_MSI.out.tsv,
                VEP.out.tsv,
                VEP.out.vcf,
                VEP.out.report,
                CALLSET_QC_TSV.out.tsv,
                // Metrics
                PICARD.out.specific_metrics,
                PICARD.out.algn_metrics,
                MOSDEPTH.out.all,
                BCFTOOLS_STATS.out.plots,
            ],
        )
        .map { meta, report, normal_bam, tumor_bam, normal_trimmed, tumor_trimmed, normal_fastp, tumor_fastp, mutect2_unfiltered, mutect2_filtered, mutect2_stats, mutect2_contamination, mutect2_read_orientation, mutect2_pileup, mutect2_segmentation, strelka2_variants, strelka2_dir, muse2_variants, muse2_txt, octopus_variants, octopus_all, clairs_variants, manta_variants, manta_indels, manta_somatic, manta_dir, delly_variants, msisensorpro_all, purecn_call, purecn_loess, purecn_cov, cnvkit_all, sigprofilerassignment_all, cnv_xref_tsv, msi_xref_tsv, vep_tsv, vep_vcf, vep_report, callset_qc_tsv, picard_specific_metrics, picard_algn_metrics, mosdepth_all, bcftools_stats_plots ->
            [
                id: meta.id,
                report: report,
                normal_bam: normal_bam,
                tumor_bam: tumor_bam,
                normal_fastp: normal_fastp,
                tumor_fastp: tumor_fastp,
                normal_trimmed: normal_trimmed,
                tumor_trimmed: tumor_trimmed,
                // Small variant callers
                mutect2_unfiltered: mutect2_unfiltered,
                mutect2_filtered: mutect2_filtered,
                mutect2_stats: mutect2_stats,
                mutect2_contamination: mutect2_contamination,
                mutect2_read_orientation: mutect2_read_orientation,
                mutect2_pileup: mutect2_pileup,
                mutect2_segmentation: mutect2_segmentation,
                strelka2_variants: strelka2_variants,
                strelka2_dir: strelka2_dir,
                muse2_variants: muse2_variants,
                muse2_txt: muse2_txt,
                octopus_variants: octopus_variants,
                octopus_all: octopus_all,
                clairs_variants: clairs_variants,
                // Other callers
                manta_variants: manta_variants,
                manta_indels: manta_indels,
                manta_somatic: manta_somatic,
                manta_dir: manta_dir,
                delly_variants: delly_variants,
                msisensorpro_all: msisensorpro_all,
                purecn_call: purecn_call,
                purecn_loess: purecn_loess,
                purecn_cov: purecn_cov,
                cnvkit_all: cnvkit_all,
                sigprofilerassignment_all: sigprofilerassignment_all,
                // Annotated calls
                cnv_xref_tsv: cnv_xref_tsv,
                msi_xref_tsv: msi_xref_tsv,
                vep_tsv: vep_tsv,
                vep_vcf: vep_vcf,
                vep_report: vep_report,
                callset_qc_tsv: callset_qc_tsv,
                // Metrics
                picard_specific_metrics: picard_specific_metrics,
                picard_algn_metrics: picard_algn_metrics,
                mosdepth_all: mosdepth_all,
                bcftools_stats_plots: bcftools_stats_plots,
            ]
        }

    to_publish_logs = Utl
        .joinFirst(
            REPORT.out.logfile,
            [
                // Small variant callers
                MUTECT2_COMPLETE.out.call_log,
                MUTECT2_COMPLETE.out.filter_log,
                MUTECT2_COMPLETE.out.contamination_log,
                MUTECT2_COMPLETE.out.read_orientation_log,
                STRELKA2.out.logfile,
                MUSE2.out.logfile,
                OCTOPUS.out.logfile,
                CLAIRS.out.logfile,
                // Other callers
                MANTA.out.logfile,
                DELLY_SV.out.logfile,
                MSISENSORPRO.out.logfile,
                PURECN_COMPLETE.out.logfile,
                CNVKIT.out.logfile,
                SIGPROFILERASSIGNMENT.out.logfile,
                // Annotated calls
                CROSS_REFERENCE_CNV.out.logfile,
                CROSS_REFERENCE_MSI.out.logfile,
                VEP.out.logfile,
                CALLSET_QC_TSV.out.logfile,
                // Metrics
                PICARD.out.logfile,
                MOSDEPTH.out.logfile,
                BCFTOOLS_STATS.out.logfile,
            ],
        )
        .map { meta, report_log, mutect2_call_log, mutect2_filter_log, mutect2_contamination_log, mutect2_read_orientation_log, strelka2_log, muse2_log, octopus_log, clairs_log, manta_log, delly_log, msisensorpro_log, purecn_log, cnvkit_log, sigprofilerassignment_log, cnv_xref_log, msi_xref_log, vep_log, callset_qc_tsv_log, picard_log, mosdepth_log, bcftools_stats_log ->
            [
                id: meta.id,
                report_log: report_log,
                // Small variant callers
                mutect2_call_log: mutect2_call_log,
                mutect2_filter_log: mutect2_filter_log,
                mutect2_contamination_log: mutect2_contamination_log,
                mutect2_read_orientation_log: mutect2_read_orientation_log,
                strelka2_log: strelka2_log,
                muse2_log: muse2_log,
                octopus_log: octopus_log,
                clairs_log: clairs_log,
                // Other callers
                manta_log: manta_log,
                delly_log: delly_log,
                msisensorpro_log: msisensorpro_log,
                purecn_log: purecn_log,
                cnvkit_log: cnvkit_log,
                sigprofilerassignment_log: sigprofilerassignment_log,
                // Annotated calls
                cnv_xref_log: cnv_xref_log,
                msi_xref_log: msi_xref_log,
                vep_log: vep_log,
                callset_qc_tsv_log: callset_qc_tsv_log,
                // Metrics
                picard_log: picard_log,
                mosdepth_log: mosdepth_log,
                bcftools_stats_log: bcftools_stats_log,
            ]
        }

    publish:
    merged_results = to_publish
    logs = to_publish_logs
    multiqc = MULTIQC.out.html.mix(MULTIQC.out.data)
    purecn_baits = PURECN_COMPLETE.out.bait_intervals
}


output {
    multiqc {
    }
    purecn_baits {
        path { "purecn_ref/" }
    }

    merged_results {
        path { s ->
            // Reports
            s.report >> "${s.id}/"

            // Inputs
            s.normal_fastp >> "${s.id}/0-preprocessed/"
            s.tumor_fastp >> "${s.id}/0-preprocessed/"
            s.normal_trimmed >> "${s.id}/0-preprocessed/"
            s.tumor_trimmed >> "${s.id}/0-preprocessed/"
            s.normal_bam >> "${s.id}/1-alignments/"
            s.tumor_bam >> "${s.id}/1-alignments/"

            // Small variant callers
            s.mutect2_unfiltered >> "${s.id}/2-variant_calls/mutect2/"
            s.mutect2_filtered >> "${s.id}/2-variant_calls/mutect2/"
            s.mutect2_stats >> "${s.id}/2-variant_calls/mutect2/"
            s.mutect2_contamination >> "${s.id}/2-variant_calls/mutect2/"
            s.mutect2_read_orientation >> "${s.id}/2-variant_calls/mutect2/"
            s.mutect2_pileup >> "${s.id}/2-variant_calls/mutect2/"
            s.mutect2_segmentation >> "${s.id}/2-variant_calls/mutect2/"
            s.strelka2_variants >> "${s.id}/2-variant_calls/strelka2/"
            s.strelka2_dir >> "${s.id}/2-variant_calls/strelka2/"
            s.muse2_variants >> "${s.id}/2-variant_calls/muse2/"
            s.muse2_txt >> "${s.id}/2-variant_calls/muse2/"
            s.octopus_variants >> "${s.id}/2-variant_calls/octopus/"
            s.octopus_all >> "${s.id}/2-variant_calls/octopus/"
            s.clairs_variants >> "${s.id}/2-variant_calls/clairs/"
            s.sigprofilerassignment_all >> "${s.id}/2-variant_calls/sigprofiler/"

            // Other calls
            s.manta_variants >> "${s.id}/2-sv_calls/manta/"
            s.manta_indels >> "${s.id}/2-sv_calls/manta/"
            s.manta_somatic >> "${s.id}/2-sv_calls/manta/"
            s.manta_dir >> "${s.id}/2-sv_calls/manta/"
            s.delly_variants >> "${s.id}/2-sv_calls/delly/"
            s.msisensorpro_all >> "${s.id}/2-sv_calls/msisensorpro/"
            s.purecn_call >> "${s.id}/2-sv_calls/purecn/"
            s.purecn_loess >> "${s.id}/2-sv_calls/purecn/"
            s.purecn_cov >> "${s.id}/2-sv_calls/purecn/"
            s.cnvkit_all >> "${s.id}/2-sv_calls/cnvkit/"

            // Annotated calls
            s.cnv_xref_tsv >> "${s.id}/3-annotations/"
            s.msi_xref_tsv >> "${s.id}/3-annotations/"
            s.vep_tsv >> "${s.id}/3-annotations/"
            s.vep_vcf >> "${s.id}/3-annotations/"
            s.vep_report >> "${s.id}/3-annotations/"
            s.callset_qc_tsv >> "${s.id}/3-annotations/"

            // Metrics
            s.picard_specific_metrics >> "${s.id}/4-metrics/"
            s.picard_algn_metrics >> "${s.id}/4-metrics/"
            s.mosdepth_all >> "${s.id}/4-metrics/"
            s.bcftools_stats_plots >> "${s.id}/4-metrics/"
        }

        index {
            path 'index.csv'
            header true
        }
    }

    logs {
        path { s ->
            s.report_log >> "${s.id}/logs/"

            // Small variant callers
            s.mutect2_call_log >> "${s.id}/logs/"
            s.mutect2_filter_log >> "${s.id}/logs/"
            s.mutect2_contamination_log >> "${s.id}/logs/"
            s.mutect2_read_orientation_log >> "${s.id}/logs/"
            s.strelka2_log >> "${s.id}/logs/"
            s.muse2_log >> "${s.id}/logs/"
            s.octopus_log >> "${s.id}/logs/"
            s.clairs_log >> "${s.id}/logs/"

            // Other callers
            s.manta_log >> "${s.id}/logs/"
            s.delly_log >> "${s.id}/logs/"
            s.msisensorpro_log >> "${s.id}/logs/"
            s.purecn_log >> "${s.id}/logs/"
            s.cnvkit_log >> "${s.id}/logs/"
            s.sigprofilerassignment_log >> "${s.id}/logs/"

            // Annotated calls
            s.cnv_xref_log >> "${s.id}/logs/"
            s.msi_xref_log >> "${s.id}/logs/"
            s.vep_log >> "${s.id}/logs/"
            s.callset_qc_tsv_log >> "${s.id}/logs/"

            // Metrics
            s.picard_log >> "${s.id}/logs/"
            s.mosdepth_log >> "${s.id}/logs/"
            s.bcftools_stats_log >> "${s.id}/logs/"
        }
    }
}
