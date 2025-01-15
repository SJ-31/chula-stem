include { whole_exome } from "./workflows/wes.nf"
include { whole_exome_tumor_only } from "./workflows/wes_tumor_only.nf"
include { rnaseq } from "./workflows/rnaseq.nf"
include { panel_of_normals } from "./workflows/panel_of_normals.nf"

workflow {
    if (params.routine == "wes" && params.tumor_only) {
        whole_exome_tumor_only()
    } else if (params.routine == "wes") {
        whole_exome()
    } else if (params.routine == "pon") {
        // Data type can be "wes" or "wgs"
        panel_of_normals()
    } else if (params.routine == "rnaseq") {
        rnaseq()
    } else {
        println "No routine specified, exiting..."
    }
}
