include { MUTECT2_COMPLETE } from "../subworkflows/mutect2_complete.nf"
include { PREPROCESS_FASTQ } from "../subworkflows/preprocess_fastq.nf"
include { PURECN_COMPLETE } from "../subworkflows/purecn_complete.nf"
include { MANTA } from "../modules/manta.nf"
include { STRELKA2 } from "../modules/strelka2.nf"
include { MSISENSORPRO } from "../modules/msisensorpro.nf"
include { MSISENSORPRO_COLLECT } from "../modules/msisensorpro.nf"
include { CNVKIT } from "../modules/cnvkit.nf"
include { CNVKIT_PREP } from "../modules/cnvkit_prep.nf"
include { MOSDEPTH } from "../modules/mosdepth.nf"
include { REPORT } from "../modules/report.nf"
include { BCFTOOLS_STATS } from "../modules/bcftools_stats.nf"
include { MULTIQC } from "../modules/multiqc.nf"
include { DELLY_SV } from "../modules/delly_sv.nf"
include { PICARD } from "../modules/picard.nf"
include { CLASSIFY_CNV_FORMAT } from "../modules/classify_cnv_format.nf"
include { OCTOPUS } from "../modules/octopus.nf"
include { CLAIRS } from "../modules/clairs.nf"
include { CLASSIFY_CNV } from "../modules/classify_cnv.nf"
include { CONCAT_VCF as CONCAT_SMALL_1 } from "../modules/concat_vcf.nf"
include { CONCAT_VCF as CONCAT_SMALL_2 } from "../modules/concat_vcf.nf"
include { CONCAT_VCF as CONCAT_SV } from "../modules/concat_vcf.nf"
include { CALLSET_QC as QC_SMALL } from "../modules/callset_qc.nf"
include { CALLSET_QC_TSV } from "../modules/callset_qc_tsv.nf"
include { CROSS_REFERENCE as CROSS_REFERENCE_MSI } from "../modules/cross_reference.nf"
include { CROSS_REFERENCE as CROSS_REFERENCE_CNV } from "../modules/cross_reference.nf"
include { STANDARDIZE_VCF } from "../modules/standardize_vcf.nf"
include { MUSE2 } from "../modules/muse2.nf"
include { GRIDSS } from "../modules/gridss.nf"
include { SIGPROFILERASSIGNMENT } from "../modules/sigprofilerassignment.nf"
include { SIGPROFILERASSIGNMENT_COLLECT } from "../modules/sigprofilerassignment.nf"
include { VEP } from "../modules/vep.nf"
include { DELLY_CNV } from '../modules/delly_cnv.nf'
include { PANEL_OF_NORMALS_FROM_BAM } from "../subworkflows/panel_of_normals_from_bam.nf"


workflow whole_exome {

    main:
    def cohort_name = params.cohort ?: "cohort"
    def branchSources = branchCriteria { it ->
        tumor: it[0].type == "cancer" || it[0].type == "tumor"
        normal: it[0].type == "normal"
    }
    def toConcat = { suffix, outdir_name, it ->
        [it[0] + ["suffix": suffix,
                  "out": "${params.outdir}/${it[0].id}/${outdir_name}",
                  "log": "${params.logdir}/${it[0].id}/${outdir_name}"]] + [it[1..-1]]
    }
    def nullIfNotNum = { it -> it.text.isNumber() ? it.text : null }
    def newOutPath = { ch, path ->
        [ch[0] + ["out": "${params.outdir}/${ch[0].id}/${path}",
                   "log": "${params.logdir}/${ch[0].id}/${path}"]] + ch[1..-1] }

    def cohortTopLevel = { it -> [["out": params.outdir, "log": params.logdir,
                             "filename": cohort_name], it] }
    /*
     * Preprocessing
     */

    PREPROCESS_FASTQ(params.input, params.outdir, params.logdir, "wes", 0)
    // After preprocessing, branch data into tumor and normal then pair up by id and join with  index


    
    branched = PREPROCESS_FASTQ.out.bam.branch(branchSources)
    normals = Utl.getId(branched.normal
                        .map({ it -> [["type": "paired",
                                 "id": it[0].id,
                                 "RGSM_normal": it[0].RGSM,
                                 "filename": it[0].id],
                                it[1]] })
                        .map({ it -> newOutPath(it, "paired")}), true)
    tumors = branched.tumor.map({ it -> [it[0].id, ["RGSM_tumor": it[0].RGSM], it[1]] })
    indices = Utl.getId(PREPROCESS_FASTQ.out.bam_index).groupTuple()
    paired = normals.join(tumors)
        .map({ it -> [it[0]] + [it[1] + it[3]] + [it[2]] + [it[4]] }) // Merge the maps
        .join(indices)
    // Order is important for the first two
    paired_no_id = Utl.delId(paired)
    // paired_no_id is [meta, normal, tumor, [normal_index, tumor_index]]
    // the order of indices doesn't matter, they just need to be present

    /*
     * (Optional) Panel of normals
     */
    if (params.create_pon) {
        to_pon_bam = branched.normal
            .map({ it -> [it[0] + ["out": "${params.outdir}/for_panel_of_normals"]] +
                  it[1..-1] }) 
        to_pon_bam_indices = PREPROCESS_FASTQ.out.bam_index
            .filter({ v -> v[0].type == "normal" })
        
        to_pon_bam.dump(tag: "pon_bam")
        to_pon_bam_indices.dump(tag: "pon_bam_indices")
        PANEL_OF_NORMALS_FROM_BAM(to_pon_bam,
                                  to_pon_bam_indices,
                                  cohort_name,
                                  channel.empty(),
                                  false)
        panel_of_normals = PANEL_OF_NORMALS_FROM_BAM.out.pon_ws.map({ it -> it[1] })
    } else {
        panel_of_normals = params.ref.panel_of_normals ?: ""
    }
    
    /*
     * Variant calling
     */

    // Structural variants

    MANTA(paired_no_id, params.ref.genome, params.ref.targets, 5)
    DELLY_SV(paired_no_id, params.ref.genome, params.ref.delly_exclude, 5)
    MSISENSORPRO(paired_no_id, params.ref.homopolymers_microsatellites, "exome",
                 params.ref.genome_gff, 5)

    MSISENSORPRO_COLLECT(MSISENSORPRO.out.summary.collect().map(cohortTopLevel), 6)
    // GRIDSS(paired_no_id, params.ref.genome, params.ref.genome_blacklist, 5)

    // Small variants
    to_mutect = paired_no_id.map({
        it -> [it[0] + ["out": "${it[0].out}/5-Mutect2"]] + it[1..-1]
    })
    MUTECT2_COMPLETE(to_mutect, panel_of_normals, 5)
    to_strelka = Utl.delId(paired.join(MANTA.out.indels))
    STRELKA2(to_strelka, params.ref.genome, params.ref.targets, 5)
    MUSE2(paired_no_id, params.ref.genome, params.ref.dbsnp, "exome", 5)

    // Combine variants by type
    small_variants_to_geno = Utl.joinFirst(MUTECT2_COMPLETE.out.filtered,
                                           [STRELKA2.out.variants
                                            .map({ it -> [it[0]] + it[1] }),
                                            MUSE2.out.variants])
        .map({ it -> toConcat("Concat_to_oct", "paired", it) })
    CONCAT_SMALL_1(small_variants_to_geno, 6)

    // Octopus and Clairs uses previous variants to aid calling
    to_geno = Utl.joinFirst(paired_no_id, [CONCAT_SMALL_1.out.vcf])
    OCTOPUS(to_geno, params.ref.genome, params.ref.targets, false, 5)
    CLAIRS(to_geno, params.ref.genome, params.ref.targets, 5)

    to_concat_small_2 = Utl.joinFirst(CONCAT_SMALL_1.out.vcf,
                                      [OCTOPUS.out.variants,
                                       CLAIRS.out.variants])
        .map({ it -> [it[0]] + [it[1..-1]]})
    CONCAT_SMALL_2(to_concat_small_2, 6)
    structural_variants = Utl.joinFirst(MANTA.out.somatic,
                                        [DELLY_SV.out.variants,
                                         // GRIDSS.out.variants
        ])
        .map({ it -> toConcat("SV_all", "annotations", it) })
    CONCAT_SV(structural_variants, 6)
    to_std_small = Utl.joinFirst(CONCAT_SMALL_2.out.vcf,
                                 [branched.normal, branched.tumor])

    /*
     * QC, Re-compute AD, DP and VAF for small variants (some callers do not compute this innately)
     */

    STANDARDIZE_VCF(Utl.addSuffix(to_std_small, "Small_std"), params.ref.genome, 6)

    small_all = Utl.delSuffix(STANDARDIZE_VCF.out.vcf)
    sv_all = Utl.delSuffix(CONCAT_SV.out.vcf)

    QC_SMALL(Utl.addSuffix(small_all, "Small_high_conf"), params.small_qc,
             panel_of_normals, 7)

    /*
    * Copy number aberration
    */

    if (!params.ref.cnvkit_reference) {
        collected_normals = normals.map({ it -> it[2] })
            .mix(PREPROCESS_FASTQ.out.bam_index.map({ it -> it[1] })).toList()
        CNVKIT_PREP(channel.of(["filename": cohort_name,
                                "out": "${params.outdir}/cnvkit_reference",
                                "log": "${params.outdir}/cnvkit_reference"])
                    .merge(collected_normals) { meta, bams -> tuple(meta, tuple(bams)) },
                    params.ref.genome, params.ref.baits_unzipped,
                    params.ref.genome_blacklist, true, "hybrid", 4)
        cnvkit_reference = CNVKIT_PREP.out.reference.first()
        cnvkit_autobin = CNVKIT_PREP.out.autobin.first()
    } else {
        cnvkit_reference = params.ref.cnvkit_reference
        cnvkit_autobin = params.ref.cnvkit_autobin
    }

    branched_with_index = PREPROCESS_FASTQ.out.bam_with_index.branch(branchSources)
    PURECN_COMPLETE(MUTECT2_COMPLETE.out.unfiltered,
                    branched_with_index.tumor,
                    branched_with_index.normal,
                    panel_of_normals,
                    cohort_name,
                    4)
    
    to_cnvkit = Utl.delId(paired.map({ it -> it[0..1] + [it[3]] })
                          .join(Utl.getId(MUTECT2_COMPLETE.out.filtered))
                          .join(Utl.getId(PURECN_COMPLETE.out.purity_ploidy)))

    CNVKIT(to_cnvkit, cnvkit_reference, "hybrid", 5)


    CLASSIFY_CNV_FORMAT(CNVKIT.out.cns, 5)
    cnv_bed = CLASSIFY_CNV_FORMAT.out.bed
        .collectFile( { meta, file -> [ "5-${meta.id}-ClassifyCNV_all.bed", file ] },
                     keepHeader: true, skip: 1)
        .map({ it -> def id = (it.baseName =~ /[0-9]+-(.*)-ClassifyCNV.*/)[0][1]
              [["filename": id, "id": id], it]
            }).map({ it -> newOutPath(it, "annotations") })

    /*
     * Variant annotation
     */
    to_vep_small = Utl.modifyMeta(small_all,
                                  ["suffix": "VEP_small",
                                   "variant_class": "small",
                                   "qc": params.small_qc])
        .map({ it -> newOutPath(it, "annotations") })
    to_vep_sv = Utl.modifyMeta(sv_all,
                               ["suffix": "VEP_SV", "variant_class": "sv",
                                "qc": params.sv_qc])

    VEP(to_vep_small.mix(to_vep_sv), params.ref.genome, panel_of_normals, 7)

    to_qc_tsv = Utl.joinFirst(VEP.out.tsv,
                              [MSISENSORPRO.out.tsv.mix(MSISENSORPRO.out.tsv)])

    CALLSET_QC_TSV(to_qc_tsv, "", 8)

    SIGPROFILERASSIGNMENT(Utl.delSuffix(QC_SMALL.out.vcf), true,
                          "${params.configdir}/excluded_signatures.txt", 7)
    SIGPROFILERASSIGNMENT_COLLECT(SIGPROFILERASSIGNMENT.out.activities
                                  .map({ it -> it[1] })
                                  .collect().map(cohortTopLevel), 8)
    CLASSIFY_CNV(cnv_bed, 7)

    // Cross reference regions
    CROSS_REFERENCE_CNV(Utl.addSuffix(CLASSIFY_CNV.out.tsv, "CR_CNV"),
                        "CNV", params.ref.clingen_dosage,
                        params.ref.cnv_reference, 8)

    to_cr_msi = Utl.modifyMeta(MSISENSORPRO.out.tsv, [suffix: "CR_MSI"])
        .map({ it -> newOutPath(it, "annotations") })
    CROSS_REFERENCE_MSI(to_cr_msi, "MSI", params.ref.clingen_gene,
                        params.ref.msi_reference, 8)

    /*
     * Metric collection
     */

    to_metrics = PREPROCESS_FASTQ.out.bam_with_index
        .map({ it -> newOutPath(it, "metrics") })

    PICARD(to_metrics, "hs", params.ref.genome,
           params.ref.targets_il, params.ref.baits_il,
           "", "", 8)
    MOSDEPTH(to_metrics, params.ref.targets, 8)

    to_bcftools_stats = Utl.addSuffix(small_all, "Bcftools_stats_small")
        .mix(Utl.addSuffix(sv_all, "Bcftools_stats_SV"))

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
                   civic_cache: params.ref.civic_cache ?: "civic_cache.json",
                   pandrugs2_cache: params.ref.pandrugs2_cache ?: "pandrugs2_cache.json",
                   cosmic_reference: params.ref.cosmic_reference]]})

    to_report = Utl.delId(Utl.getId(others, true).join(vep_to_report))
        .map({ it -> [it[0], it[1] + it[2]] }).map({ newOutPath(it, "") })

    caches = []
    REPORT(to_report, caches, "variant_calling", 8)

}
