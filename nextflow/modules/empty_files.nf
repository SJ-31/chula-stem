process EMPTY_FILES {
    ext version: "1"
    // Generate empty files to meet process input requirements when they don't actually need it

    publishDir "${meta.out}", mode:"copy", saveAs: params.saveFn
    publishDir "${meta.log}", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(file)
    val(count)
    //

    output:
    tuple val(meta), path("*.empty")
    //

    script:
    name = file.baseName
    """
    for i in \$(seq 1 ${count}); do
        touch "${name}_\${i}.empty"
    done
    """
    //
}

// To unpack the paths i.e. [meta, [f1, f2, f3, ...]] -> [meta, f1, f2, f3, ...]
// Use the following closure:
// EMPTY_FILES.out.map({ [it[0] + it[1]] })
