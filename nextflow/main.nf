include { whole_exome } from "./subworkflows/wes"

workflow {
    switch(params.routine) {
        case "wes":
            whole_exome();
            break;
    }
}
