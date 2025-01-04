include { PREPROCESS_FASTQ } from "../subworkflows/preprocess_fastq.nf"
include { MUTECT2_COMPLETE } from "../subworkflows/mutect2_complete.nf"
include { SIGPROFILERASSIGNMENT } from "../modules/sigprofilerassignment.nf"
include { MANTA } from "../modules/manta.nf"
include { MSISENSORPRO } from "../modules/msisensorpro.nf"
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
include { CLAIRS } from '../modules/clairs.nf'
include { MUTECT2 } from '../modules/mutect2.nf'

workflow whole_exome_tumor_only {

    main:
    def cohort_name = params.cohort ? params.cohort : "cohort"
    /*
     * Preprocessing
     */
    PREPROCESS_FASTQ(params.input, params.outdir, params.logdir, "wes")

    empty_normals = EMPTY_FILES_1(PREPROCESS_FASTQ.out.bam, 1)
    tumors = PREPROCESS_FASTQ.out.bam.map { [
                        ["type": "paired",
                        "id": it[0].id,
                        "out": "${params.outdir}/${it[0].id}/variant_calling",
                        "filename": it[0].id,
                         "RGSM_tumor": it[0].RGSM,
                        "log": "${params.logdir}/${it[0].id}/variant_calling"],
                                        it[1]] }

    indices = PREPROCESS_FASTQ.out.bam_index
    empty_indices = EMPTY_FILES_2(PREPROCESS_FASTQ.out.bam_index, 1)
    all_indices = indices.map(Utl.getId).mix(empty_indices.map(Utl.getId))
        .groupTuple()

    paired = Utl.joinFirst(empty_normals, [tumors])
        .map({ [it[0]] + [it[2]] + [it[1]] + [it[3]] })
        .join(all_indices)

    paired_no_id = paired.map(Utl.delId)

    /*
     * Variant calling
     */
    // Structural variants

    MANTA(paired_no_id, params.ref.genome, params.ref.targets, 5)
    MSISENSORPRO(paired_no_id, params.ref.homopolymers_microsatellites, "exome",
                 params.ref.genome_gff, 5)
    GRIDSS(paired_no_id, params.ref.genome, params.ref.genome_blacklist, 5)

    // Small variants

    to_mutect = paired_no_id.map { [it[0] + ["out": "${it[0].out}/5-Mutect2"]] + it[1..-1] }
    MUTECT2_COMPLETE(to_mutect, 5)

    to_clairs = Utl.joinFirst(tumors, [PREPROCESS_FASTQ.out.bam_index,
                                       MUTECT2_COMPLETE.out])

    CLAIRS_TO(to_clairs, params.ref.genome, params.ref.targets, 5)

    def toConcat = { suffix, outdir_name, it ->
        [it[0] + ["suffix": suffix,
                  "out": "${params.outdir}/${it[0].id}/${outdir_name}",
                  "log": "${params.logdir}/${it[0].id}/${outdir_name}"]] + [it[1..-1]]
    }

    // Combine variants by type
    small_variants_to_oct = Utl.joinFirst(
        MUTECT2_COMPLETE.out, [CLAIRS_TO.out.variants]
    ).map({ toConcat("Small_all", "annotations", it) })

    CONCAT_SMALL_1(small_variants_to_oct, params.ref.genome, 6)

    to_octopus = Utl.joinFirst(paired_no_id, [CONCAT_SMALL_1.out.vcf])

    OCTOPUS(to_octopus, params.ref.genome, params.ref.targets, 5)

    CONCAT_SMALL_2(Utl.joinFirst(CONCAT_SMALL_1.out.vcf, [OCTOPUS.out.variants])
                    .map({ toConcat("Small_all", "annotations", it) }),
                   params.ref.genome, 6)

    small_variants = CONCAT_SMALL_2.out.vcf

    structural_variants = Utl.joinFirst(MANTA.out.somatic, [GRIDSS.out.variants])
        .map({ toConcat("SV_all", "annotations", it) })

    CONCAT_SV(structural_variants, params.ref.genome, 6)

    /*
     * Metric collection and QC
     */

    to_standardize = Utl.addSuffix(Utl.joinFirst(small_variants,
                                                 [empty_normals,
                                                  PREPROCESS_FASTQ.out.bam]),
                                   "Small_std")

    STANDARDIZE_VCF(to_standardize, params.ref.genome, 6)

    QC_SMALL(Utl.addSuffix(STANDARDIZE_VCF.out.vcf,
                           "Small_high_conf"), params.small_qc, 7)

    QC_SV(Utl.addSuffix(CONCAT_SV.out.vcf, "SV_high_conf"), params.sv_qc, 7)

    /*
    * Copy number abberation
    */
    purity_ploidy = tumors.map({ [it[0].id, null, null] })

    CNVKIT_PREP(Channel.of(["filename": "flat_reference",
                            "out": "${params.outdir}/cnvkit_cnn",
                            "log": "${params.outdir}/cnvkit_cnn"])
                    .merge(Channel.fromPath("${params.configdir}/EMPTY.txt")),
                params.ref.genome, params.ref.baits_unzipped,
                params.ref.genome_blacklist, 4)

    to_cnvkit = paired.map({ it[0..1] + [it[3]] })
            .join(QC_SMALL.out.vcf.map(Utl.getId))
            .join(purity_ploidy)
            .map(Utl.delId)

    CNVKIT(to_cnvkit, CNVKIT_PREP.out.reference.first(), "hybrid", 5)

    cnv_ch = CNVKIT.out.cnv

    CLASSIFY_CNV_FORMAT(cnv_ch, 5)
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
    // QC for VEP will be carried out separately
    to_vep_small = Utl.modifyMeta(STANDARDIZE_VCF.out.vcf,
                              ["suffix": "VEP_small", "variant_class": "small",
                               "qc": params.small_qc])


    to_vep_sv = Utl.modifyMeta(CONCAT_SV.out.vcf,
                           ["suffix": "VEP_SV", "variant_class": "sv",
                            "qc": params.sv_qc])


    VEP(to_vep_small.mix(to_vep_sv), params.ref.genome, 7)

    vep_out = VEP.out.tsv.branch { meta, files ->
        sv: meta.variant_class == "sv"
        small: meta.variant_class == "small"
    }

    to_qc_tsv_sv = Utl.joinFirst(vep_out.sv, [MSISENSORPRO.out.tsv])
    to_qc_tsv_small = Utl.joinFirst(vep_out.small, [MSISENSORPRO.out.tsv])

    CALLSET_QC_TSV(to_qc_tsv_sv.mix(to_qc_tsv_small), "", 8)

    SIGPROFILERASSIGNMENT(Utl.delSuffix(QC_SMALL.out.vcf), true,
                          "${params.configdir}/excluded_signatures.txt", 7)
    CLASSIFY_CNV(cnv_bed, 7)

    CROSS_REFERENCE_CNV(CLASSIFY_CNV.out.tsv, "CNV", params.ref.clingen_dosage,
                        params.ref.cnv_reference, 8)
    CROSS_REFERENCE_MSI(MSISENSORPRO.out.tsv, "MSI", params.ref.clingen_gene,
                        params.ref.msi_reference, 8)

    /*
     * Metric collection
     */
    def replaceOut = { [it[0] + ["out": "${params.outdir}/${it[0].id}/metrics",
                                 "log": "${params.logdir}/${it[0].id}/metrics"],
                        it[1], it[2]] }

    // def getIdType = {  [it[0].id, it[0].type, it[0], it[1]]  }

    // to_metrics = PREPROCESS_FASTQ.out.bam.map(getIdType)
    //     .join(PREPROCESS_FASTQ.out.bam_index.map(getIdType),
    //           by: [0, 1]) // [id, type, meta, bam, meta, bai]
    //     .map({ [it[2], it[3], it[5]] })
    //     .map(replaceOut)
    //     TODO: <2025-01-04 Sat> Check if this works
    to_metrics = Utl.joinFirst(PREPROCESS_FASTQ.out.bam,
                               [PREPROCESS_FASTQ.out.bam_index],
                               on: ["id", "type"]).map(replaceOut)
    // <2025-01-04 Sat> If joinFirst is set up correctly, then "id" "type" that
    // was prepended should be removed automatically

    PICARD(to_metrics, "hs", params.ref.genome,
           params.ref.targets_il, params.ref.baits_il,
           "", 8)
    MOSDEPTH(to_metrics, params.ref.targets, 8)

    to_bcftools_stats = Utl.delSuffix(STANDARDIZE_VCF.out.vcf)
        .mix(Utl.delSuffix(CONCAT_SV.out.vcf))

    BCFTOOLS_STATS(to_bcftools_stats, params.ref.targets, 8)

    to_multiqc = PREPROCESS_FASTQ.out.fastp_json.mix(VEP.out.report,
                                                    MOSDEPTH.out.dist,
                                                    PICARD.out.metrics,
                                                    BCFTOOLS_STATS.out.py)
        .flatten().collect().map({ [["out": params.outdir, "log": params.logdir,
        "filename": cohort_name], it] })
    MULTIQC(to_multiqc, 8)

}

