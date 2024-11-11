process VEP {
    ext version: ""
    conda { task.ext.conda }

    publishDir "$meta.out", mode: "copy"
    publishDir "$meta.log", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(vcf)
    //

    output:
    path("*.log")
    //

    script:
    output =
        check = file("meta.out/output")
    if (check.exists()) {
        """

        cp meta.log/ .
        """
    } else {
        """


        cp .command.out .log
        """
    }
    //
}
