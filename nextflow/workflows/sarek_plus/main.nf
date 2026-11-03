include { OPTITYPE } from "./modules/optitype.nf"
include { PURECN_COMPLETE } from "./subworkflows/purecn_complete.nf"

workflow {

    main:
    input = channel.fromPath(workflow.outputDir)
    // TODO: build the manifest file from the outputDir, which is the
    // same output directory of sarek.
    // Throw an error if there aren't enough samples for it to work
    // TODO: use optitype to get the hla classification given fastp

    
    OPTITYPE(input.map { val -> [val, [val.fastp_1, val.fastp2]] }, "dna")
    mutect2_unfiltered = input.map { val -> [val, val.mutect2_unfiltered] }
    
    PURECN_COMPLETE(mutect2_unfiltered)
}
