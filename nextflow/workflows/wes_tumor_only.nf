include { PREPROCESS_FASTQ } from "../subworkflows/preprocess_fastq.nf"
include { MUTECT2_COMPLETE } from "../subworkflows/mutect2_complete.nf"
include { SIGPROFILERASSIGNMENT } from "../modules/sigprofilerassignment.nf"
include { SIGPROFILERASSIGNMENT_COLLECT } from "../modules/sigprofilerassignment.nf"
include { MANTA } from "../modules/manta.nf"
include { MSISENSORPRO } from "../modules/msisensorpro.nf"
include { MSISENSORPRO_COLLECT } from "../modules/msisensorpro.nf"
include { EMPTY_FILES as EMPTY_FILES_1 } from "../modules/empty_files.nf"
include { EMPTY_FILES as EMPTY_FILES_2 } from "../modules/empty_files.nf"
include { CNVKIT } from "../modules/cnvkit.nf"
include { CNVKIT_PREP } from "../modules/cnvkit_prep.nf"
include { MOSDEPTH } from "../modules/mosdepth.nf"
include { BCFTOOLS_STATS } from "../modules/bcftools_stats.nf"
include { MULTIQC } from "../modules/multiqc.nf"
include { PICARD } from "../modules/picard.nf"
include { CONCAT_VCF as CONCAT_SMALL_1 } from "../modules/concat_vcf.nf"
include { CONCAT_VCF as CONCAT_SMALL_2 } from "../modules/concat_vcf.nf"
include { CONCAT_VCF as CONCAT_SV } from "../modules/concat_vcf.nf"
include { CLAIRS_TO } from "../modules/clairs_to.nf"
include { CLASSIFY_CNV_FORMAT } from "../modules/classify_cnv_format.nf"
include { CLASSIFY_CNV } from "../modules/classify_cnv.nf"
include { CALLSET_QC_TSV } from "../modules/callset_qc_tsv.nf"
include { OCTOPUS } from "../modules/octopus.nf"
include { CALLSET_QC as QC_SMALL } from "../modules/callset_qc.nf"
include { CALLSET_QC as QC_SV } from "../modules/callset_qc.nf"
include { STANDARDIZE_VCF } from "../modules/standardize_vcf.nf"
include { GRIDSS } from "../modules/gridss.nf"
include { CROSS_REFERENCE as CROSS_REFERENCE_MSI } from "../modules/cross_reference.nf"
include { CROSS_REFERENCE as CROSS_REFERENCE_CNV } from "../modules/cross_reference.nf"
include { FACETS_PILEUP } from "../modules/facets_pileup.nf"
include { FACETS } from "../modules/facets.nf"
include { VEP } from "../modules/vep.nf"
include { CLAIRS } from "../modules/clairs.nf"
include { MUTECT2 } from "../modules/mutect2.nf"
include { REPORT } from "../modules/report.nf"
include { PURECN_COMPLETE } from "../subworkflows/purecn_complete.nf" 
include { GET_THERAPY_CACHE } from "../modules/get_therapy_cache.nf"

workflow whole_exome_tumor_only {

    main:
    def cohort_name = params.cohort ? params.cohort : "cohort"
    def cohortTopLevel = { it -> [["out": params.outdir, "log": params.logdir,
                             "filename": cohort_name], it] }

    def panel_of_normals = params.ref.panel_of_normals ? params.ref.panel_of_normals : ""
    /*
     * Preprocessing
     */
    PREPROCESS_FASTQ(params.input, params.outdir, params.logdir, "wes", 0)

    empty_normals = EMPTY_FILES_1(PREPROCESS_FASTQ.out.bam, "", 1)
    tumors = PREPROCESS_FASTQ.out.bam.map { it -> [
                        ["type": "paired",
                        "id": it[0].id,
                        "out": "${params.outdir}/${it[0].id}/variant_calling",
                        "filename": it[0].id,
                         "RGSM_tumor": it[0].RGSM,
                        "log": "${params.logdir}/${it[0].id}/variant_calling"],
                                        it[1]] }
    // tumors: [meta, bam]

    indices = PREPROCESS_FASTQ.out.bam_index
    empty_indices = EMPTY_FILES_2(PREPROCESS_FASTQ.out.bam_index, "", 1)
    all_indices = Utl.getId(indices).mix(Utl.getId(empty_indices))
        .groupTuple()

    paired = Utl.joinFirst(empty_normals, [tumors],
                           ["id"], true) // -> [n_meta, n, t_meta, t]
        .map({ it -> [it[2].id] + [it[2]] + [it[1]] + [it[3]] }) // -> [id, t_meta, n, t]
        .join(all_indices) // -> [id, t_meta, n, t, [n_index, t_index]]

    paired_no_id = Utl.delId(paired)

    /*
     * Variant calling
     */
    // Structural variants

    MANTA(paired_no_id, params.ref.genome, params.ref.targets, 5)
    MSISENSORPRO(paired_no_id, params.ref.homopolymers_microsatellites, "exome",
                 params.ref.genome_gff, 5)
    MSISENSORPRO_COLLECT(MSISENSORPRO.out.summary.collect().map(cohortTopLevel), 6)
    GRIDSS(paired_no_id, params.ref.genome, params.ref.genome_blacklist, 5)

    // Small variants

    to_mutect = paired_no_id.map { it -> [it[0] + ["out": "${it[0].out}/5-Mutect2"]] + it[1..-1] }
    MUTECT2_COMPLETE(to_mutect, panel_of_normals, 5)

    to_clairs = Utl.joinFirst(tumors, [PREPROCESS_FASTQ.out.bam_index,
                                       MUTECT2_COMPLETE.out.filtered])

    CLAIRS_TO(to_clairs, params.ref.genome, params.ref.targets, 5)

    def toConcat = { suffix, outdir_name, it ->
        [it[0] + ["suffix": suffix,
                  "out": "${params.outdir}/${it[0].id}/${outdir_name}",
                  "log": "${params.logdir}/${it[0].id}/${outdir_name}"]] + [it[1..-1]]
    }

    // Combine variants by type
    small_variants_to_oct = Utl.joinFirst(
        MUTECT2_COMPLETE.out.filtered, [CLAIRS_TO.out.variants]
    ).map({ it -> toConcat("Small_all", "annotations", it) })

    CONCAT_SMALL_1(small_variants_to_oct, 6)

    to_octopus = Utl.joinFirst(paired_no_id, [CONCAT_SMALL_1.out.vcf])

    OCTOPUS(to_octopus, params.ref.genome, params.ref.targets, true, 5)

    CONCAT_SMALL_2(Utl.joinFirst(CONCAT_SMALL_1.out.vcf, [OCTOPUS.out.variants])
                    .map({ it -> toConcat("Small_all", "annotations", it) }), 6)

    small_variants = CONCAT_SMALL_2.out.vcf

    structural_variants = Utl.joinFirst(MANTA.out.somatic, [GRIDSS.out.variants])
        .map({ it -> toConcat("SV_all", "annotations", it) })

    CONCAT_SV(structural_variants, 6)

    /*
     * Metric collection and QC
     */

    to_standardize = Utl.addSuffix(Utl.joinFirst(small_variants,
                                                 [empty_normals,
                                                  PREPROCESS_FASTQ.out.bam]),
                                   "Small_std")

    STANDARDIZE_VCF(to_standardize, params.ref.genome, 6)

    QC_SMALL(Utl.addSuffix(STANDARDIZE_VCF.out.vcf, "Small_high_conf"),
             params.small_qc,
             panel_of_normals,
             7)

    QC_SV(Utl.addSuffix(CONCAT_SV.out.vcf, "SV_high_conf"), params.sv_qc,
          panel_of_normals, 7)

    /*
    * Copy number abberation
     */
    if (params.ref.purecn_normaldb) {
        PURECN_COMPLETE(MUTECT2_COMPLETE.out.unfiltered,
                        PREPROCESS_FASTQ.out.bam_with_index,
                        channel.empty(),
                        4)
        purity_ploidy = Utl.getId(PURECN_COMPLETE.out.purity_ploidy)
    } else {
        purity_ploidy = tumors.map({ it -> [it[0].id, null, null] })
    }
    
    if (!params.ref.cnvkit_reference) {
        CNVKIT_PREP(channel.of(["filename": "flat_reference",
                                "out": "${params.outdir}/cnvkit_cnn",
                                "log": "${params.outdir}/cnvkit_cnn"])
                        .merge(channel.fromPath("${params.configdir}/EMPTY.txt")),
                    params.ref.genome, params.ref.baits_unzipped,
                    params.ref.genome_blacklist, false, "hybrid", 4)
        cnvkit_reference = CNVKIT_PREP.out.reference.first()
    } else {
        cnvkit_reference = params.ref.cnvkit_reference
    }

    to_cnvkit = Utl.delId(paired.map({ it -> it[0..1] + [it[3]] })
                          .join(Utl.getId(QC_SMALL.out.vcf))
                          .join(purity_ploidy))

    CNVKIT(to_cnvkit, cnvkit_reference, "hybrid", 5)

    CLASSIFY_CNV_FORMAT(CNVKIT.out.cns, 5)

    cnv_bed = CLASSIFY_CNV_FORMAT.out.bed
        .collectFile( { it -> [ "5-${it[0].id}-ClassifyCNV_all.bed", it[1] ] },
                     keepHeader: true, skip: 1)
        .map({ it -> def id = (it.baseName =~ /[0-9]+-(.*)-ClassifyCNV.*/)[0][1]
              [["id": id, "filename": id,
          "out": "${params.outdir}/${id}/annotations",
          "log": "${params.logdir}/${id}/annotations"], it]
        })
    /*
     * Variant annotation
     */
    // QC for VEP will be carried out separately on the TSV file,
    // so it does not receive the QC vcf
    to_vep_small = Utl.modifyMeta(STANDARDIZE_VCF.out.vcf,
                              ["suffix": "VEP_small", "variant_class": "small",
                               "qc": params.small_qc])

    to_vep_sv = Utl.modifyMeta(CONCAT_SV.out.vcf,
                           ["suffix": "VEP_SV", "variant_class": "sv",
                            "qc": params.sv_qc])

    VEP(to_vep_small.mix(to_vep_sv), params.ref.genome, panel_of_normals, 7)

    to_qc_tsv = Utl.joinFirst(VEP.out.tsv,
                              [MSISENSORPRO.out.tsv.mix(MSISENSORPRO.out.tsv)])

    CALLSET_QC_TSV(to_qc_tsv, "", 8)

    SIGPROFILERASSIGNMENT(Utl.delSuffix(QC_SMALL.out.vcf), true,
                          "${params.configdir}/excluded_signatures.txt", 7)
    SIGPROFILERASSIGNMENT_COLLECT(SIGPROFILERASSIGNMENT.out.activities.map({ it -> it[1] })
                                    .collect().map(cohortTopLevel), 8)
    CLASSIFY_CNV(cnv_bed, 7)

    // Cross reference regions
    CROSS_REFERENCE_CNV(Utl.addSuffix(CLASSIFY_CNV.out.tsv, "CR_CNV"),
                        "CNV", params.ref.clingen_dosage,
                        params.ref.cnv_reference, 8)
    to_cr_msi = Utl.modifyMeta(MSISENSORPRO.out.tsv,
                               [suffix: "CR_MSI",
                                out: { it -> "${params.outdir}/${it.id}/annotations" },
                                log: { it -> "${params.logdir}/${it.id}/annotations" } ])
    CROSS_REFERENCE_MSI(to_cr_msi, "MSI", params.ref.clingen_gene,
                        params.ref.msi_reference, 8)

    /*
     * Metric collection
     */
    def replaceOut = { it -> [it[0] + ["out": "${params.outdir}/${it[0].id}/metrics",
                                 "log": "${params.logdir}/${it[0].id}/metrics"],
                        it[1], it[2]] }

    to_metrics = Utl.joinFirst(PREPROCESS_FASTQ.out.bam,
                               [PREPROCESS_FASTQ.out.bam_index],
                               ["id", "type"]).map(replaceOut)

    PICARD(to_metrics, "hs", params.ref.genome,
           params.ref.targets_il, params.ref.baits_il,
           "", "", 8)
    MOSDEPTH(to_metrics, params.ref.targets, 8)


    to_bcftools_stats = Utl.addSuffix(STANDARDIZE_VCF.out.vcf, "Bcftools_stats_small")
        .mix(Utl.addSuffix(CONCAT_SV.out.vcf, "Bcftools_stats_SV"))

    BCFTOOLS_STATS(to_bcftools_stats, params.ref.targets, 8)

    to_multiqc = PREPROCESS_FASTQ.out.fastp_json.mix(VEP.out.report,
                                                    MOSDEPTH.out.dist,
                                                    PICARD.out.metrics,
                                                    BCFTOOLS_STATS.out.py)
        .flatten().collect().map(cohortTopLevel)
    MULTIQC(to_multiqc, 8)

    // Combine channels for report
    vep_grouped = Utl.getId(CALLSET_QC_TSV.out.tsv)
        .groupTuple(sort: { it -> (it.baseName =~ /VEP_small/) ? 1 : -1 })
    vep_to_report = vep_grouped
        .map({ it -> [it[0]] + [vep_small: it[1][1], vep_sv: it[1][0]] })

    others = Utl.joinFirst(SIGPROFILERASSIGNMENT.out.activities,
                           [CNVKIT.out.cns, CNVKIT.out.cnr,
                            CROSS_REFERENCE_CNV.out.tsv,
                            CROSS_REFERENCE_MSI.out.tsv ]).map(
        { it -> [it[0], [sigprofiler: it[1], cnvkit_cns: it[2],
                   cnvkit_cnr: it[3], classify_cnv: it[4],
                   msisensor_pro: it[5],
                   civic_cache: params.ref.civic_cache ? params.ref.civic_cache : "civic_cache.json",
                   pandrugs2_cache: params.ref.pandrugs2_cache ? params.ref.pandrugs2_cache : "pandrugs2_cache.json",
                   cosmic_reference: params.ref.cosmic_reference], ]})

    to_report = Utl.modifyMeta(Utl.delId(Utl.getId(others, true).join(vep_to_report))
                               .map({ it -> [it[0], it[1] + it[2]] }),
                               [out: { it -> "${params.outdir}/${it.id}" },
                                log: { it -> "${params.logdir}/${it.id}" }])

    // if (!params.ref.civic_cache && !params.ref.pandrugs2_cache) {
    //     all_vep = vep_grouped.map({ it -> it[1] }).flatten().collect().map({
    //         it -> [[out: "${params.outdir}/cache", log: params.logdir], it]
    //     })
    //     GET_THERAPY_CACHE(all_vep)
    //     caches = GET_THERAPY_CACHE.out.civic
    //         .mix(GET_THERAPY_CACHE.out.pandrugs2).collect()
    // } else {
    //     caches = channel.fromPath(params.ref.civic_cache)
    //         .mix(channel.fromPath(params.ref.pandrugs2_cache)).collect()
    // }
    caches = [] // <2025-01-08 Wed> for temporary debugging
    REPORT(to_report, caches, "variant_calling", 8)
}
