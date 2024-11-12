include { whole_exome } from "./workflows/wes"

workflow {
    switch(params.routine) {
        case "wes":
            whole_exome();
            break;
    }
}
