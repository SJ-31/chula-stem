process PURECN_BAIT_INTERVALS {
    ext version: ""

    publishDir "${meta.out}", mode:"copy", saveAs: params.saveFn
    publishDir "${meta.log}", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), val(genome)
    val(baits)
    val(module_number)
    //

    output:
    path("*.log")
    //

    script:
    output = Utl.getName()
    check = file("${meta.out}/${output}")
    args = task.ext.args.join(" ")
    if (check.exists()) {
        """
        ln -sr ${check} .
        ln -sr ${meta.log}/ .
        """
    } else {
        """
        

        cp .command.out .log
        """
    }
    //
}
