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
include { FACETS_PILEUP } from "../modules/facets_pileup.nf"
include { FACETS } from "../modules/facets.nf"
include { VEP } from "../modules/vep.nf"


workflow whole_exome_tumor_only {

    main:
    def cohort_name = params.cohort ? params.cohort : "cohort"
    /*
     * Preprocessing
     */
    PREPROCESS_FASTQ(params.input, params.outdir, params.logdir)

    empty_normals = EMPTY_FILES_1(PREPROCESS_FASTQ.out.bam, 1)
    tumors = PREPROCESS_FASTQ.out.bam.map { [it[0].id,
                        ["type": "paired",
                        "id": it[0].id,
                        "out": "${params.outdir}/${it[0].id}/variant_calling",
                        "filename": it[0].id,
                         "RGSM_tumor": it[0].RGSM,
                        "log": "${params.logdir}/${it[0].id}/variant_calling"],
                                        it[1]] }

    indices = PREPROCESS_FASTQ.out.bam_index
    empty_indices = EMPTY_FILES_2(PREPROCESS_FASTQ.out.bam_index, 1)
    all_indices = indices.map(params.getId).mix(empty_indices.map(params.getId))
        .groupTuple()

    paired = empty_normals.map(params.getId).join(tumors)
        .map({ [it[0]] + [it[2]] + [it[1]] + [it[3]] })
        .join(all_indices)

    paired_no_id = paired.map(params.delId)

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

    to_clairs = tumors.join(PREPROCESS_FASTQ.out.bam_index.map(params.getId))
            .join(MUTECT2_COMPLETE.out.map(params.getId))
            .map(params.delId)
    CLAIRS_TO(to_clairs, params.ref.genome, params.ref.targets, 5)

    def toConcat = { suffix, outdir_name, it ->
        [it[0] + ["suffix": suffix,
                  "out": "${params.outdir}/${it[0].id}/${outdir_name}",
                  "log": "${params.logdir}/${it[0].id}/${outdir_name}"]] + [it[1..-1]]
    }

    // Combine variants by type
    small_variants_to_oct = MUTECT2_COMPLETE.out.map(params.prependId)
        .join(CLAIRS_TO.out.variants.map(params.getId))
        .map(params.delId)
        .map({ toConcat("Small_all", "annotations", it) })
        .join(PREPROCESS_FASTQ.out.bam.map(params.getId))

    CONCAT_SMALL_1(small_variants_to_oct, params.ref.genome, 6)

    to_octopus = paired_no_id.join(CONCAT_SMALL_1.out.vcf.map(params.getId))
    OCTOPUS(to_octopus, params.ref.genome, params.ref.targets, 5)

     CONCAT_SMALL_2(CONCAT_SMALL_1.out.vcf.map(params.prependId)
                    .join(OCTOPUS.out.variants.map(params.getId)),
                    params.ref.genome, 6)
    small_variants = CONCAT_SMALL_2.out.vcf

    structural_variants = MANTA.out.somatic.map(params.prependId)
        .join(GRIDSS.out.variants.map(params.getId))
        .map(params.delId)
        .map({ toConcat("SV_all", "annotations", it) })

    CONCAT_SV(structural_variants, params.ref.genome, 6)

    /*
     * Metric collection and QC
     */
    STANDARDIZE_VCF(small_variants.map({ params.addSuffix("Small_std", it) }),
                    params.ref.genome, 6)

    QC_SMALL(STANDARDIZE_VCF.out.vcf.map({ params.addSuffix("Small_high_conf", it) }),
        params.small_qc, 7)
    QC_SV(CONCAT_SV.out.vcf.map({ params.addSuffix("SV_high_conf", it) }),
          params.sv_qc, 7)

    /*
    * Copy number abberation
    */

    // FACETS_PILEUP(paired_no_id, params.ref.pileup, 5) // TODO: this is only temporary,
    // // in the real run use dbSNP or the combined snps file instead
    // FACETS(FACETS_PILEUP.out.pileup, 5)

    // TODO: only need this with facets
    // def nullIfNotNum = { it.text.isNumber() ? it.text : null }

    purity_ploidy = paired.map({ [it[0], null, null] })
    // purity_ploidy = FACETS.out.purity_ploidy
    //     .map({ [it[0].id, nullIfNotNum(it[1]), nullIfNotNum(it[2])] })

    CNVKIT_PREP(Channel.of(["filename": "flat_reference",
                            "out": "${params.outdir}/cnvkit_cnn",
                            "log": "${params.outdir}/cnvkit_cnn"])
                    .merge(Channel.fromPath("${params.configdir}/EMPTY.txt")),
                params.ref.genome, params.ref.baits_unzipped,
                params.ref.genome_blacklist, 4)

    to_cnvkit = paired.map({ it[0..1] + [it[2]] })
            .join(QC_SMALL.out.vcf.map(params.getId))
            .join(purity_ploidy)
            .map(params.delId)

    CNVKIT(to_cnvkit, CNVKIT_PREP.out.reference, "hybrid", 5)

    cnv_ch = CNVKIT.out.cnv
    // TODO: maybe add in facets?

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
    vep_small_add = ["suffix": "VEP_small", "variant_class": "small", "qc": params.small_qc]
    to_vep_small = STANDARDIZE_VCF.out.vcf.map({ params.addMeta(vep_small_add, it) })

    vep_sv_add = ["suffix": "VEP_SV", "variant_class": "sv", "qc": params.sv_qc]
    to_vep_sv = CONCAT_SV.out.vcf.map({ params.addMeta(vep_sv_add, it) })

    VEP(to_vep_small.mix(to_vep_sv), params.ref.genome, 7)
    def joinTwo = { x, y ->
        x.map(params.prependId).join(y.map(params.getId)).map(params.delId)
    }

    vep_out = VEP.out.tsv.branch { meta, files ->
        sv: meta.variant_class == "sv"
        small: meta.variant_class == "small"
    }
    to_qc_tsv_sv = joinTwo(vep_out.sv, MSISENSORPRO.out.tsv)
    to_qc_tsv_small = joinTwo(vep_out.small, MSISENSORPRO.out.tsv)
    CALLSET_QC_TSV(to_qc_tsv_sv.mix(to_qc_tsv_small), "", 8)

    SIGPROFILERASSIGNMENT(QC_SMALL.out.vcf.map({ params.addSuffix(null, it) }), true,
                          "${params.configdir}/excluded_signatures.txt", 7)
    CLASSIFY_CNV(cnv_bed, 7)

    /*
     * Metric collection
     */
    def replaceOut = { [it[0] + ["out": "${params.outdir}/${it[0].id}/metrics",
                                 "log": "${params.logdir}/${it[0].id}/metrics"],
                        it[1], it[2]] }

    def getIdType = {  [it[0].id, it[0].type, it[0], it[1]]  }
    to_metrics = PREPROCESS_FASTQ.out.bam.map(getIdType)
        .join(PREPROCESS_FASTQ.out.bam_index.map(getIdType),
              by: [0, 1]) // [id, type, meta, bam, meta, bai]
        .map({ [it[2], it[3], it[5]] })
        .map(replaceOut)

    PICARD(to_metrics, "hs", params.ref.genome, params.ref.targets_il, params.ref.baits_il,
           "", 8)
    MOSDEPTH(to_metrics, params.ref.targets, 8)
    BCFTOOLS_STATS(STANDARDIZE_VCF.out.vcf.mix(CONCAT_SV.out.vcf), params.ref.targets, 8)

    to_multiqc = PREPROCESS_FASTQ.out.fastp_json.mix(VEP.out.report,
                                                    MOSDEPTH.out.dist,
                                                    PICARD.out.metrics,
                                                    BCFTOOLS_STATS.out.py)
        .flatten().collect().map({ [["out": params.outdir, "log": params.logdir,
        "filename": cohort_name], it] })
    MULTIQC(to_multiqc, 8)

}
