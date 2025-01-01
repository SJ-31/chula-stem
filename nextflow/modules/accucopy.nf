process ACCUCOPY {
    ext version: "c58cf6bd-7L0B04WU-debug"

    publishDir "$meta.out", mode:"copy", saveAs: params.saveFn
    publishDir "$meta.log", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(normal), path(tumor), path(indices, arity: "2")
    val(reference)
    val(variants) // Common snps associated with the reference, in BED format
    val(module_number)
    //

    output:
    path("*.log")
    //

    script:
    output = "${module_number}-${meta.filename}-Accucopy"
    check = file("${meta.out}/${output}")
    split = reference.split("\\.")
    basename = split.size() > 1 ? (split - split[-1]).join(".") : split[0]
    args = task.ext.args.join(" ")
    def configure_text = """read_length\t${task.ext.read_length}
    window_size\t${task.ext.window_size}
    reference_folder_path\treference
    samtools_path\t/usr/local/bin/samtools
    caller_path\t/usr/local/strelka
    accucopy_path\t/usr/local/Accucopy"""
    if (check.exists()) {
        """
        cp -r ${check}
        ln -sr ${meta.log}/accucopy.log .
        """
    } else {
        """
        mkdir reference
        ln -sr ${reference} reference/genome.fa
        ln -sr "${reference}.fai" reference/genome.fa.fai
        ln -sr "${basename}.dict" reference/genome.dict
        ln -sr ${variants} reference/snp_sites.gz
        ln -sr "${variants}.tbi" reference/snp_sites.gz.tbi

        echo -e "${configure_text}" > configure.tsv

        /usr/local/Accucopy/main.py \\
            --configure_filepath configure.tsv \\
            --tumor_bam ${tumor} \\
            --normal_bam ${normal} \\
            --output_dir ${output}

        get_nextflow_log.bash accucopy.log
        """
    }
    //
}
