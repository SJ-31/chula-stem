include { OPTITYPE } from "./modules/optitype.nf"
include { MERGE_FASTQ } from "./modules/merge_fastq.nf"
// include { PURECN_COMPLETE } from "./subworkflows/purecn_complete.nf"
include { SAREK_MANIFEST } from "./modules/sarek_manifest.nf"

nextflow.enable.types = true
params {
    sarek_config: Path
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

    merged.view()

    fastq = preprocessed.filter { v -> !v.split_fastq }.map { v -> record(id: v.id, fastq_1: file(v.fastq_1), fastq_2: file(v.fastq_2), status: v.status) }.mix(merged)

    optitype = OPTITYPE(fastq.filter { v -> v.status == 0 }, "dna")

    publish:

    // TODO: use 'topics' to extract versions as well as logs
    // https://docs.seqera.io/nextflow/tutorials/topic-channels/

    // TODO: extract VEP annotations from the VCF
    //      sarek uses "CSQ" as the VEP field instead of "ANN"
    //

    // TODO: need a routine "extract depth" that retrieves allele depths
    // and other statistics for every caller
    // [2026-11-05 Thu] Written this already, just need to add it now

    // TODO: [2026-11-05 Thu] after parsing the output directory
    // and manifest, will produce three channels: germline, paired, and
    // tumor-only channels
    // The latter go to PURECN and structural variant calling
    // The rest go straight to the report
    // `paired` has columns for tumor and normal recal files

    // OPTITYPE(germline.map { val -> [val, [val.fastp_1, val.fastp2]] }, "dna")

    // mutect2_unfiltered = input.map { val -> [val, val.mutect2_unfiltered] }

    // // TODO: add a step that combines structural variant callers and
    // passes them to vep. Just use concat_vcf

    // PURECN_COMPLETE(mutect2_unfiltered)
    optitype = optitype
}

output {
    optitype {
        path { v ->
            v.optitype_csv >> "optitype/${v.id}"
            v.optitype_pdf >> "optitype/${v.id}"
        }
    }
}
