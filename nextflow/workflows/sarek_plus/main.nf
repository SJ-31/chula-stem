include { OPTITYPE } from "./modules/optitype.nf"
include { MERGE_FASTQ } from "./modules/merge_fastq.nf"
include { PURECN_COMPLETE } from "./subworkflows/purecn_complete.nf"
include { SAREK_MANIFEST } from "./modules/sarek_manifest.nf"
include { CONCAT_VCF } from "./modules/concat_vcf.nf"
include { SIGPROFILERASSIGNMENT ; SIGPROFILERASSIGNMENT_COLLECT } from "./modules/sigprofilerassignment.nf"
include { EXTRACT_FORMAT_FIELDS } from "./modules/extract_format_fields.nf"
include { FORMAT_VEP_VCF } from "./modules/format_vep_vcf.nf"

nextflow.enable.types = true
params {
    sarek_config: Path
    ref: Map
    defaults: Map
}

workflow {

    main:
    sarek_out = SAREK_MANIFEST(params.sarek_config)

    preprocessed = sarek_out.preprocessed
        .flatMap { csv ->
            csv.splitCsv(header: true)
        }
        .map { row ->
            record(id: row.id, status: row.status as Integer, recal: row.recal, recal_index: row.recal_index, fastq_1: row.fastq_1, fastq_2: row.fastq_2, split_fastq: row.split_fastq)
        }

    merged = MERGE_FASTQ(preprocessed.filter { v -> v.split_fastq }.map { v -> record(id: v.id, fastq_dir: v.fastq_1, status: v.status) })

    fastq = preprocessed.filter { v -> !v.split_fastq }.map { v -> record(id: v.id, fastq_1: file(v.fastq_1), fastq_2: file(v.fastq_2), status: v.status) }.mix(merged)

    fastq.filter { v -> v.status == 0 }.view()

    optitype = OPTITYPE(fastq.filter { v -> v.status == 0 }, "dna")

    to_purecn = sarek_out.tool_outs
        .flatMap { csv ->
            csv.splitCsv(header: true)
        }
        .map { row ->
            record(id: row.id, mutect2_vcf: row.mutect2_unfiltered)
        }
        .join(preprocessed.filter { v -> v.status == 0 }.map { v -> record(id: v.id, normal_bam: v.recal, normal_bam_index: v.recal_index) }, by: "id")
        .join(preprocessed.filter { v -> v.status == 1 }.map { v -> record(id: v.id, tumor_bam: v.recal, tumor_bam_index: v.recal_index) }, by: "id")
    to_purecn.view()

    sv_tumor_to_concat = sarek_out.tool_outs.filter { v -> v.output_type == "paired" || v.output_type == "tumor_only" }.map { v -> record(id: v.id, vcfs: [v?.vep_manta_somatic, v?.vep_tiddit]) }

    sv_tumor_concat = CONCAT_VCF(sv_tumor_to_concat, "sv_")

    to_extract_format = sarek_out.tool_outs.map {
        { v ->
            def res = [:]
            ["deepvariant", "freebayes", "strelka", "tiddit", "deepsomatic", "manta_somatic", "muse", "mutect2", "strelka_indels", "haplotypecaller", "deepvariant"].each { el ->
                if (v[el]) {
                    res[el] = v[el]
                }
            }
            record(id: v.id, caller2vcfs: res, suffix: "_${v.output_type}", normal_sample: v.normal_sample, tumor_sample: v.tumor_sample)
        }
    }
    to_extract_format.view()
    format_fields = EXTRACT_FORMAT_FIELDS(to_extract_format)

    to_vep_tsv = sarek_out.tool_outs.filter { v -> v.vep_consensus }.map { v -> record(id: v.id, vcf: v.vep_consensus) }.mix(sv_tumor_concat.map { v -> record(id: v.id, vcf: v.concat) })
    vep_tsvs = FORMAT_VEP_VCF(to_vep_tsv, ["CALLERS", "NCALLERS", "SVTYPE"], "CSQ")

    publish:

    // TODO: use 'topics' to extract versions as well as logs
    // https://docs.seqera.io/nextflow/tutorials/topic-channels/

    // TODO: extract VEP annotations from the VCF
    //      sarek uses "CSQ" as the VEP field instead of "ANN"
    //

    // TODO: need to merge the output of the CNV callers

    optitype = optitype
    format_fields = format_fields
    vep_tsvs = vep_tsvs
}

output {
    optitype {
        path { v ->
            v.optitype_csv >> "optitype/${v.id}"
            v.optitype_pdf >> "optitype/${v.id}"
        }
    }
    vep_tsvs {
        path { v -> v.tsv >> "annotation/${v.id}" }
    }
    format_fields {
        path { v -> v.tsv >> "format_fields/${v.id}" }
    }
}
