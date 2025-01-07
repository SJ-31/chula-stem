include { whole_exome } from "./workflows/wes.nf"
include { whole_exome_tumor_only } from "./workflows/wes_tumor_only.nf"

workflow {
    if (params.routine == "wes" && params.tumor_only) {
        whole_exome_tumor_only()
    } else if (params.routine == "wes") {
            whole_exome()
    } else {
        println "No routine specified, exiting..."
    }
}
