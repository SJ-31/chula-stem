include { FASTP } from "../modules/fastp.nf"
include { BWA } from "../modules/bwa.nf"
include { MARK_DUPLICATES } from "../modules/mark_duplicates.nf"
include { BQSR } from "../modules/bqsr.nf"
include { MUTECT2_COMPLETE } from "../subworkflows/mutect2_complete.nf"
include { MANTA } from "../modules/manta.nf"
include { STRELKA2 } from "../modules/strelka2.nf"
include { MSISENSORPRO } from "../modules/msisensorpro.nf"
include { CNVKIT } from "../modules/cnvkit.nf"
include { DELLY_SV } from "../modules/delly_sv.nf"
include { PICARD } from "../modules/picard.nf"
include { CONCAT_VCF as CONCAT_SMALL } from "../modules/concat_vcf.nf"
include { CONCAT_VCF as CONCAT_SV } from "../modules/concat_vcf.nf"
include { CALLSET_QC as QC_SMALL } from "../modules/callset_qc.nf"
include { STANDARDIZE_VCF } from "../modules/standardize_vcf.nf"
include { MUSE2 } from "../modules/muse2.nf"
include { GRIDSS } from "../modules/gridss.nf"
include { FACETS_PILEUP } from "../modules/facets_pileup.nf"
include { FACETS } from "../modules/facets.nf"
include { DELLY_COV } from "../modules/delly_cov.nf"
include { VEP } from "../modules/vep.nf"
include { SAMTOOLS_INDEX } from "../modules/samtools_index.nf"

workflow "whole_exome" {

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

    FASTP(input, 1)
    BWA(FASTP.out.passed, params.ref.genome, 2)
    MARK_DUPLICATES(BWA.out.mapped, params.ref.genome, 3)
    BQSR(MARK_DUPLICATES.out.dedup, params.ref.genome, params.ref.known_variants, 4)

    // After preprocessing, branch data into tumor and normal then pair up by id and join with
    //  index
    output = BQSR.out.bam
    SAMTOOLS_INDEX(output) // Required by certain callers
    branched = output.branch { meta, files ->
        tumor: meta.type == "cancer" || meta.type == "tumor"
        normal: meta.type == "normal"
    }
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

    def toConcat = { suffix, outdir_name, variant_class, it ->
        [it[0] + ["suffix": suffix,
                  "vclass": variant_class,
                  "out": "${params.outdir}/${it[0].id}/${outdir_name}",
                  "log": "${params.logdir}/${it[0].id}/${outdir_name}"]] + [it[1..-1]]
    }

    // Combine variants by type
    small_variants = MUTECT2_COMPLETE.out.map(params.prependId)
        .join(STRELKA2.out.variants.map(params.getId).map( { it.flatten() } ))
        .join(MUSE2.out.variants.map(params.getId))
        .map({ it[1..-1] })
        .map({ toConcat("Small_std", "annotations", "small", it) })

    CONCAT_SMALL(small_variants, 6)

    structural_variants = MANTA.out.somatic.map(params.prependId)
        .join(DELLY_SV.out.variants.map(params.getId))
        .join(GRIDSS.out.variants.map(params.getId))
        .map({ it[1..-1] })
        .map({ toConcat("SV_std", "annotations", "SV", it) })

    CONCAT_SV(structural_variants, 6)

    def withBams = { n, t, it -> it.map(params.prependId)
                    .join(n.map(params.getId))
                    .join(t.map(params.getId))
                    .map( { it[1..-1] } ) }

    std_small_variants = withBams(branched.normal, branched.tumor, CONCAT_SMALL.out.vcf)
    std_svs = withBams(branched.normal, branched.tumor, CONCAT_SV.out.vcf)
    to_standardize = std_small_variants.mix(std_svs)

    /*
     * Metric collection and QC
     */
    // TODO: concat, get basic stats, then run qc on the whole thing before effect prediction
    // CONCAT_SMALL.out.vcf
    STANDARDIZE_VCF(to_standardize, params.ref.genome, 7)

    standardized = STANDARDIZE_VCF.out.vcf.branch { meta, files ->
        small: meta.vclass == "small"
        SV: meta.vclass == "SV"
    }

    // TODO: qc must be done differently for SVs

    // QC_SMALL(standardized.small.map({ addSuffix("Small_high_conf", it) }),
        // params.small_qc, 8)
    // callset_qc(???, ???, ???, ???)
    //

    /*
    * Copy number abberation
    */

    FACETS_PILEUP(paired_no_id, params.ref.pileup, 5) // TODO: this is only temporary,
    // in the real run use dbSNP or the combined snps file instead
    FACETS(FACETS_PILEUP.out.pileup, 5)

    def nullIfNotNum = { it.text.isNumber() ? it.text : null }
    purity_ploidy = FACETS.out.purity_ploidy
        .map({ [it[0], nullIfNotNum(it[1]), nullIfNotNum(it[2])] })

    to_delly_cnv = paired.join(DELLY_COV.out.cov).join(purity_ploidy).map({it[1..-1]})
    // DELLY_CNV(to_delly_cnv, params.ref.genome, params.ref.delly_mappability, 5)

    // collected_normals = normals.map({ it[1] }).collect()
    // CNVKIT_PREP([["filename": "cohort"], collected_normals], params.ref.baits,
    //             params.ref.genome_blacklist, 4)

    // final_small = ??? TODO: the small variants after callset qc
    // to_cnvkit = paired.join(final_small, purity_ploidy)
    // CNVKIT(to_cnvkit, CNVKIT_PREP.out.reference, "wgs", 5)

    /*
     * Variant annotation
     */
    // VEP(???, params.ref.genome, 7)
    //
    // TODO: need to transfer annotations between them
    //      Check if merge1 works for this
    // annotated = SNPEFF.out.vcf.join(VEP.map(getId))
    // MERGE_B(annotated, 8)


    // PICARD("???", "exome", params.ref.genome, params.ref.targets, params.ref.baits,
    //      null, 8)

}
