include { MUTECT2_COMPLETE } from "../subworkflows/mutect2_complete.nf"
include { PREPROCESS_FASTQ } from "../subworkflows/preprocess_fastq.nf"
include { MANTA } from "../modules/manta.nf"
include { STRELKA2 } from "../modules/strelka2.nf"
include { MSISENSORPRO } from "../modules/msisensorpro.nf"
include { CNVKIT } from "../modules/cnvkit.nf"
include { CNVKIT_PREP } from "../modules/cnvkit_prep.nf"
include { MOSDEPTH } from "../modules/mosdepth.nf"
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
// include { CALLSET_QC as QC_SV } from "../modules/callset_qc.nf"
include { CALLSET_QC_TSV } from "../modules/callset_qc_tsv.nf"
include { CROSS_REFERENCE as CROSS_REFERENCE_MSI } from "../modules/cross_reference.nf"
include { CROSS_REFERENCE as CROSS_REFERENCE_CNV } from "../modules/cross_reference.nf"
include { STANDARDIZE_VCF } from "../modules/standardize_vcf.nf"
include { MUSE2 } from "../modules/muse2.nf"
include { GRIDSS } from "../modules/gridss.nf"
include { FACETS_PILEUP } from "../modules/facets_pileup.nf"
include { FACETS } from "../modules/facets.nf"
include { SIGPROFILERASSIGNMENT } from "../modules/sigprofilerassignment.nf"
include { VEP } from "../modules/vep.nf"
include { DELLY_CNV } from '../modules/delly_cnv.nf'


workflow whole_exome {

    main:
    def cohort_name = params.cohort ? params.cohort : "cohort"
    def branchSources = branchCriteria { meta, files ->
        tumor: meta.type == "cancer" || meta.type == "tumor"
        normal: meta.type == "normal"
    }
    /*
     * Preprocessing
     */
    PREPROCESS_FASTQ(params.input, params.outdir, params.logdir, "wes", 0)
    // After preprocessing, branch data into tumor and normal then pair up by id and join with
    //  index
    branched = PREPROCESS_FASTQ.out.bam.branch(branchSources)
    normals = branched.normal.map { [it[0].id,
                                        ["type": "paired",
                                        "id": it[0].id,
                                        "out": "${params.outdir}/${it[0].id}/paired",
                                        "RGSM_normal": it[0].RGSM,
                                        "filename": it[0].id,
                                        "log": "${params.logdir}/${it[0].id}/paired"],
                                        it[1]] }
    tumors = branched.tumor.map { [it[0].id, ["RGSM_tumor": it[0].RGSM], it[1]] }
    indices = Utl.getId(PREPROCESS_FASTQ.out.bam_index).groupTuple()
    paired = normals.join(tumors)
        .map({ [it[0]] + [it[1] + it[3]] + [it[2]] + [it[4]] }) // Merge the maps
        .join(indices)
    // Order is important for the first two
    paired_no_id = Utl.delId(paired)
    // paired_no_id is [meta, normal, tumor, [normal_index, tumor_index]]
    // the order of indices doesn't matter, they just need to be present

    /*
     * Variant calling
     */
    // Structural variants

    MANTA(paired_no_id, params.ref.genome, params.ref.targets, 5)
    DELLY_SV(paired_no_id, params.ref.genome, params.ref.delly_exclude, 5)
    MSISENSORPRO(paired_no_id, params.ref.homopolymers_microsatellites, "exome",
                 params.ref.genome_gff, 5)
    GRIDSS(paired_no_id, params.ref.genome, params.ref.genome_blacklist, 5)

    // Small variants

    to_mutect = paired_no_id.map { [it[0] + ["out": "${it[0].out}/5-Mutect2"]] + it[1..-1] }
    MUTECT2_COMPLETE(to_mutect, 5)

    to_strelka = Utl.delId(paired.join(MANTA.out.indels))

    STRELKA2(to_strelka, params.ref.genome, params.ref.targets, 5)
    MUSE2(paired_no_id, params.ref.genome, params.ref.dbsnp, "exome", 5)


    def toConcat = { suffix, outdir_name, it ->
        [it[0] + ["suffix": suffix,
                  "out": "${params.outdir}/${it[0].id}/${outdir_name}",
                  "log": "${params.logdir}/${it[0].id}/${outdir_name}"]] + [it[1..-1]]
    }


    // Combine variants by type
    //
    small_variants_to_geno = Utl.joinFirst(MUTECT2_COMPLETE.out,
                                           [STRELKA2.out.variants.map({ it.flatten }),
                                            MUSE2.out.variants])
        .map({ toConcat("Concat_to_oct", "paired", it) })

    CONCAT_SMALL_1(small_variants_to_geno, params.ref.genome, 6)

    // Octopus and Clairs uses previous variants to aid calling
    to_geno = Utl.joinFirst(paired, [CONCAT_SMALL_1.out.vcf])
    OCTOPUS(to_geno, params.ref.genome, params.ref.targets, 5)
    CLAIRS(to_geno, params.ref.genome, params.ref.targets, 5)

    to_concat_small_2 = Utl.joinFirst(CONCAT_SMALL_1.out.vcf,
                                      [OCTOPUS.out.variants,
                                       CLAIRS.out.variants])
    CONCAT_SMALL_2(to_concat_small_2, params.ref.genome, 6)

    structural_variants = Utl.joinFirst(MANTA.out.somatic,
                                        [DELLY_SV.out.variants, GRIDSS.out.variants])
        .map({ toConcat("SV_all", "annotations", it) })

    CONCAT_SV(structural_variants, params.ref.genome, 6)

    to_std_small = Utl.joinFirst(CONCAT_SMALL_2.out.vcf,
                                 [branched.normal, branched.tumor])
    /*
     * QC, Re-compute AD, DP and VAF for small variants (some callers do not compute this innately)
     */

    STANDARDIZE_VCF(Utl.addSuffix(to_std_small, "Small_std"), params.ref.genome, 6)

    small_all = Utl.delSuffix(STANDARDIZE_VCF.out.vcf)
    sv_all = Utl.delSuffix(CONCAT_SV.out.vcf)

    QC_SMALL(Utl.addSuffix(small_all, "Small_high_conf"), params.small_qc, 7)

    /*
    * Copy number abberation
    */
    FACETS_PILEUP(paired_no_id, params.ref.pileup, 5)
    // TODO: join this with an estimate of window_size from cnvkit's autobin
    FACETS(FACETS_PILEUP.out.pileup, 5)

    def nullIfNotNum = { it.text.isNumber() ? it.text : null }
    purity_ploidy = FACETS.out.purity_ploidy
        .map({ [it[0].id, nullIfNotNum(it[1]), nullIfNotNum(it[2])] })

    collected_normals = normals.map({ it[2] }).toList()
    CNVKIT_PREP(Channel.of(["filename": cohort_name,
                            "out": "${params.outdir}/cnvkit_reference",
                            "log": "${params.outdir}/cnvkit_reference"])
                            .merge(collected_normals),
                params.ref.genome, params.ref.baits_unzipped, params.ref.genome_blacklist, 4)

    to_cnvkit = Utl.delId(paired.map({ it[0..1] + [it[3]] })
            .join(Utl.getId(QC_SMALL.out.vcf))
            .join(purity_ploidy))

    CNVKIT(to_cnvkit, CNVKIT_PREP.out.reference.first(), "hybrid", 5)

    CLASSIFY_CNV_FORMAT(CNVKIT.out.cns.mix(FACETS.out.rds), 5)
    cnv_bed = CLASSIFY_CNV_FORMAT.out.bed
        .collectFile( { meta, file -> [ "5-${meta.id}-ClassifyCNV_all.bed", file ] },
                     keepHeader: true, skip: 1)
        .map({ def id = (it.baseName =~ /.*-(.*)-.*/)[0][1]
        [["filename": id,
          "out": "${params.outdir}/${id}/annotations",
          "log": "${params.logdir}/${id}/annotations"], it]
        })

    /*
     * Variant annotation
     */
    to_vep_small = Utl.modifyMeta(small_all,
                                  ["suffix": "VEP_small", "variant_class": "small",
                                   "qc": params.small_qc])
    to_vep_sv = Utl.modifyMeta(sv_all,
                               ["suffix": "VEP_SV", "variant_class": "sv",
                                "qc": params.sv_qc])

    VEP(to_vep_small.mix(to_vep_sv), params.ref.genome, 7)

    to_qc_tsv = Utl.joinFirst(VEP.out.tsv,
                              [MSISENSORPRO.out.tsv.mix(MSISENSORPRO.out.tsv)])

    CALLSET_QC_TSV(to_qc_tsv, "", 8)

    SIGPROFILERASSIGNMENT(Utl.delSuffix(QC_SMALL.out.vcf), true,
                          "${params.configdir}/excluded_signatures.txt", 7)
    CLASSIFY_CNV(cnv_bed, 7)

    // Cross reference regions
    CROSS_REFERENCE_CNV(Utl.addSuffix(CLASSIFY_CNV.out.tsv, "CR_CNV"),
                        "CNV", params.ref.clingen_dosage,
                        params.ref.cnv_reference, 8)
    to_cr_msi = Utl.modifyMeta(MSISENSORPRO.out.tsv,
                               [suffix: "CR_MSI",
                                out: { "${params.outdir}/${it.id}/annotations" },
                                log: { "${params.logdir}/${it.id}/annotations" } ])
    CROSS_REFERENCE_MSI(to_cr_msi, "MSI", params.ref.clingen_gene,
                        params.ref.msi_reference, 8)

    /*
     * Metric collection
     */
    def replaceOut = { [it[0] + ["out": "${params.outdir}/${it[0].id}/metrics",
                                 "log": "${params.logdir}/${it[0].id}/metrics"],
                        it[1], it[2]] }

    to_metrics = Utl.joinFirst(PREPROCESS_FASTQ.out.bam,
                               [PREPROCESS_FASTQ.out.bam_index],
                               on: ["id", "type"]).map(replaceOut)

    PICARD(to_metrics, "hs", params.ref.genome,
           params.ref.targets_il, params.ref.baits_il,
           "", 8)
    MOSDEPTH(to_metrics, params.ref.targets, 8)

    to_bcftools_stats = Utl.addSuffix(small_all, "Bcftools_stats_small")
        .mix(Utl.addSuffix(sv_all, "Bcftools_stats_SV"))
    BCFTOOLS_STATS(to_bcftools_stats, params.ref.targets, 8)

    to_multiqc = PREPROCESS_FASTQ.out.fastp_json.mix(VEP.out.report,
                                                    MOSDEPTH.out.dist,
                                                    PICARD.out.metrics,
                                                    BCFTOOLS_STATS.out.py)
        .flatten().collect().map({ [["out": params.outdir, "log": params.logdir,
                                     "filename": cohort_name], it] })
    MULTIQC(to_multiqc, 8)

}
