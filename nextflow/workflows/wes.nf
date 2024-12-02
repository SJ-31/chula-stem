include { FASTP } from "../modules/fastp.nf"
include { BWA } from "../modules/bwa.nf"
include { MARK_DUPLICATES } from "../modules/mark_duplicates.nf"
include { BQSR } from "../modules/bqsr.nf"
include { MUTECT2_COMPLETE } from "../subworkflows/mutect2_complete.nf"
include { MANTA } from "../modules/manta.nf"
include { STRELKA2 } from "../modules/strelka2.nf"
include { MSISENSORPRO } from "../modules/msisensorpro.nf"
include { CNVKIT } from "../modules/cnvkit.nf"
include { CNVKIT_PREP } from "../modules/cnvkit_prep.nf"
include { DELLY_SV } from "../modules/delly_sv.nf"
include { PICARD } from "../modules/picard.nf"
include { CONCAT_VCF as CONCAT_SMALL } from "../modules/concat_vcf.nf"
include { CLASSIFY_CNV_FORMAT } from "../modules/classify_cnv_format.nf"
include { CLASSIFY_CNV } from "../modules/classify_cnv.nf"
include { CONCAT_VCF as CONCAT_SV } from "../modules/concat_vcf.nf"
include { CALLSET_QC as QC_SMALL } from "../modules/callset_qc.nf"
include { CALLSET_QC as QC_SV } from "../modules/callset_qc.nf"
include { STANDARDIZE_VCF } from "../modules/standardize_vcf.nf"
include { MUSE2 } from "../modules/muse2.nf"
include { GRIDSS } from "../modules/gridss.nf"
include { FACETS_PILEUP } from "../modules/facets_pileup.nf"
include { FACETS } from "../modules/facets.nf"
include { DELLY_COV } from "../modules/delly_cov.nf"
include { VEP } from "../modules/vep.nf"
include { SAMTOOLS_INDEX } from "../modules/samtools_index.nf"


workflow whole_exome {

    main:
    input = Channel.fromPath(params.input)
        .splitCsv(header: true)
        .map { [["id": it.patient,
                "out": "${params.outdir}/${it.patient}/${it.source}",
                "type": it.source,
                "log": "${params.logdir}/${it.patient}/${it.source}",
                "filename": "${it.patient}_${it.source}", // Output filename
                "RGLB": it.read_group_library, // Read group info used by GATK tools
                "RGPL": it.read_group_platform,
                "RGPU": it.read_group_platform_unit,
                "RGSM": it.read_group_sample_name
            ], [file(it.fastq_1), file(it.fastq_2)]] }

    /*
     * Preprocessing
     */
    def branchSources = branchCriteria { meta, files ->
        tumor: meta.type == "cancer" || meta.type == "tumor"
        normal: meta.type == "normal"
    }

    FASTP(input, 1)
    BWA(FASTP.out.passed, params.ref.genome, 2)
    MARK_DUPLICATES(BWA.out.mapped, params.ref.genome, 3)
    BQSR(MARK_DUPLICATES.out.dedup, params.ref.genome, params.ref.known_variants, 4)

    // After preprocessing, branch data into tumor and normal then pair up by id and join with
    //  index
    output = BQSR.out.bam
    SAMTOOLS_INDEX(output) // Required by certain callers
    branched = output.branch(branchSources)
    normals = branched.normal.map { [it[0].id,
                                         ["type": "paired",
                                          "id": it[0].id,
                                          "out": "${params.outdir}/${it[0].id}/paired",
                                          "RGSM_normal": it[0].RGSM,
                                          "filename": it[0].id,
                                          "log": "${params.logdir}/${it[0].id}/paired"],
                                         it[1]] }
    tumors = branched.tumor.map { [it[0].id, ["RGSM_tumor": it[0].RGSM], it[1]] }
    indices = SAMTOOLS_INDEX.out.index.map({ [it[0].id, it[1]] }).groupTuple()
    paired = normals.join(tumors)
        .map({ [it[0]] + [it[1] + it[3]] + [it[2]] + [it[4]] }) // Merge the maps
        .join(indices)
    // Order is important for the first two
    paired_no_id = paired.map { it[1..-1] }
    // paired_no_id is [meta, normal, tumor, indices]

    to_delly_cov = branched.tumor.map({ [it[0].id] + it[0..-1] }).join(indices).map { it[1..-1] }
    DELLY_COV(to_delly_cov, params.ref.genome, params.ref.mappability, 4)

    /*
     * Variant calling
     */
    // Structural variants

    MANTA(paired_no_id, params.ref.genome, params.ref.targets, 5)
    DELLY_SV(paired_no_id, params.ref.genome, params.ref.delly_exclude, 5)
    MSISENSORPRO(paired_no_id, params.ref.homopolymers_microsatellites, 5)
    GRIDSS(paired_no_id, params.ref.genome, params.ref.genome_blacklist, 5)

    // Small variants

    to_mutect = paired_no_id.map { [it[0] + ["out": "${it[0].out}/5-Mutect2"]] + it[1..-1] }
    MUTECT2_COMPLETE(to_mutect, 5)

    to_strelka = paired.join(MANTA.out.indels)
        .map { it[1..-1] }
    STRELKA2(to_strelka, params.ref.genome, params.ref.targets, 5)
    MUSE2(paired_no_id, params.ref.genome, params.ref.dbsnp, "exome", 5)

    def toConcat = { suffix, outdir_name, it ->
        [it[0] + ["suffix": suffix,
                  "out": "${params.outdir}/${it[0].id}/${outdir_name}",
                  "log": "${params.logdir}/${it[0].id}/${outdir_name}"]] + [it[1..-1]]
    }

    // Combine variants by type
    small_variants = MUTECT2_COMPLETE.out.map(params.prependId)
        .join(STRELKA2.out.variants.map(params.getId).map( { it.flatten() } ))
        .join(MUSE2.out.variants.map(params.getId))
        .map({ it[1..-1] })
        .map({ toConcat("Small_all", "annotations", it) })

    CONCAT_SMALL(small_variants, params.ref.genome, 6)

    structural_variants = MANTA.out.somatic.map(params.prependId)
        .join(DELLY_SV.out.variants.map(params.getId))
        .join(GRIDSS.out.variants.map(params.getId))
        .map({ it[1..-1] })
        .map({ toConcat("SV_all", "annotations", it) })

    CONCAT_SV(structural_variants, params.ref.genome, 6)

    def withBams = { n, t, it -> it.map(params.prependId)
                    .join(n.map(params.getId))
                    .join(t.map(params.getId))
                    .map( { it[1..-1] } ) }

    std_small_variants = withBams(branched.normal, branched.tumor, CONCAT_SMALL.out.vcf)

    /*
     * Metric collection and QC
     */
    STANDARDIZE_VCF(std_small_variants.map({ params.addSuffix("Small_std", it) }),
                    params.ref.genome, 6)

    QC_SMALL(STANDARDIZE_VCF.out.vcf.map({ params.addSuffix("Small_high_conf", it) }),
        params.small_qc, 7)
    QC_SV(CONCAT_SV.out.vcf.map({ params.addSuffix("SV_high_conf", it) }),
          params.sv_qc, 7)

    /*
    * Copy number abberation
    */

    FACETS_PILEUP(paired_no_id, params.ref.pileup, 5) // TODO: this is only temporary,
    // in the real run use dbSNP or the combined snps file instead
    FACETS(FACETS_PILEUP.out.pileup, 5)

    def nullIfNotNum = { it.text.isNumber() ? it.text : null }
    purity_ploidy = FACETS.out.purity_ploidy
        .map({ [it[0].id, nullIfNotNum(it[1]), nullIfNotNum(it[2])] })

    // to_delly_cnv = paired.join(DELLY_COV.out.cov).join(purity_ploidy).map({it[1..-1]})
    // DELLY_CNV(to_delly_cnv, params.ref.genome, params.ref.delly_mappability, 5)

    collected_normals = normals.map({ it[2] }).toList()
    CNVKIT_PREP(Channel.of(["filename": "cohort",
                            "out": "${params.outdir}/cnvkit_reference",
                            "log": "${params.outdir}/cnvkit_reference"])
                            .merge(collected_normals),
                params.ref.genome, params.ref.baits_unzipped, params.ref.genome_blacklist, 4)


    to_cnvkit = paired.map({ it[0..1] + [it[2]] })
            .join(QC_SMALL.out.vcf.map(params.getId))
            .join(purity_ploidy)
            .map({ it[1..-1] })

    CNVKIT(to_cnvkit, CNVKIT_PREP.out.reference, "hybrid", 5)

    CLASSIFY_CNV_FORMAT(CNVKIT.out.cnv.mix(FACETS.out.rds), 5)
    cnv_bed = CLASSIFY_CNV_FORMAT.out.bed
        .collectFile( { meta, file -> [ "5-${meta.id}-ClassifyCNV_all.bed", file ] },
                     keepHeader: true, skip: 1)
        .map({ def id = (it.baseName =~ /.*-(.*)-.*/)[0][1]
        [["filename": id,
          "out": "${params.outdir}/${id}/annotations",
          "log": "${params.logdir}/${id}/annotations"], it]
        })
    // TODO: should have a better way of unifying CNV callers
    /*
     * Variant annotation
     */
    VEP(QC_SMALL.out.vcf.map({ params.addSuffix("VEP_small", it) })
        .mix(QC_SV.out.vcf.map({ params.addSuffix("VEP_SV", it) })),
        params.ref.genome, 7)
    CLASSIFY_CNV(cnv_bed, 7)

    /*
     * Metric collection
     */
    def replaceOut = { [it[0] + ["out": "${params.outdir}/${it[0].id}/metrics",
                                 "log": "${params.logdir}/${it[0].id}/metrics"],
                        it[1], it[2]] }

    def getIdType = {  [it[0].id, it[0].type, it[0], it[1]]  }
    to_metrics = output.map(getIdType)
        .join(SAMTOOLS_INDEX.out.index.map(getIdType), by: [0, 1]) // [id, type, meta, bam, meta, bai]
        .map({ [it[2], it[3], it[5]] })
        .map(replaceOut)
    // PICARD(to_metrics, "hs", params.ref.genome, params.ref.targets, params.ref.baits,
    //      null, 8)
    // MOSDEPTH(to_metrics, params.ref.targets, 8)
    // BCFTOOLS_STATS(STANDARDIZE_VCF.out.vcf.mix(concat_sv.out.vcf), params.ref.targets, 8)

}
