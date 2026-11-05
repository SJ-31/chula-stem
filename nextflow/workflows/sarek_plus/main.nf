include { OPTITYPE } from "./modules/optitype.nf"
include { PURECN_COMPLETE } from "./subworkflows/purecn_complete.nf"
include { SAREK_MANIFEST } from "./modules/sarek_manifest.nf"

workflow {

    main:
    input = SAREK_MANIFEST(params.sarek_config)

    germline = input.germline.splitCsv(header: true)
    tumor_only = input.tumor_only.splitCsv(header: true)
    paired = input.paired.splitCsv(header: true)

    tumor_only.view()
    
    // TODO: use optitype to get the hla classification given fastp

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
}
