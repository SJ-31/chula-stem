process SOMATICSEQ {
    ext version: ""

    publishDir "${meta.out}", mode:"copy", saveAs: params.saveFn
    publishDir "${meta.log}", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(tumor), path(normal), path(vcfs), \
        path(other_snps, arity: "1..*"), path(other_indels, arity: "1..*")
    val(reference)
    val(blacklist)
    val(target_intervals)
    val(snp_classifier)
    val(indel_classifier)
    val(module_number)
    //

    output:
    path("*.log")
    //

    script:
    def getCaller = { n, candidates ->
        def found = candidates.findAll({ (it.baseName.toLowerCase() =~ n).find() })
        if (found.size() > 1) {
            throw new Exception("Ambiguous candidates for ${n}")
        } else if (found.size() == 0) {
            throw new Exception("No candidates for ${n}")
        }
        found[0]
    }

    output = Utl.getName(module_number, "Somaticseq")
    check = file("${meta.out}/${output}")
    target_flag = target_intervals != "" ? " --inclusion-region ${target_intervals} " : ""
    blacklist_flag = blacklist != "" ? " --exclusion-region ${blacklist} " : ""

    mutect2 = getCaller("mutect2", vcfs)
    strelka_snv = getCaller("snvs_strelka", vcfs)
    strelka_indel = getCaller("indels_strelka")
    muse = getCaller("muse", vcfs)
    vcfs.remove(mutect2)
    vcfs.remove(muse)
    vcfs.remove(strelka_snv)
    vcfs.remove(strelka_indel)

    if (vcfs.size() > 0) {
        array = files.toList().collect({"\"${it}\""}).join(" ")
        remaining = (0..vcfs.size()).step(1)
        snvs = remaining.collect({ "snvs_${it}.vcf" }).join(" ")
        snvs_flag = "--arbitrary-snvs ${snvs}"
        indels = remaining.collect({ "indels_${it}.vcf" }).join(" ")
        indels_flag = "--arbitrary-indels ${indels}"
    } else {
        snvs_flag = ""
        indels_flag = ""
        array = ""
    }
    array = "(${array})"

    caller_flags = "--muse-vcf ${muse} --strelka-snv ${strelka_snv} --strelka-indel ${strelka_indel} --mutect2-vcf ${mutect2} ${indels_flag} ${snvs_flag}"
    args = task.ext.args.join(" ")
    if (check.exists()) {
        """
        ln -sr ${check} .
        ln -sr ${meta.log}/ .
        """
    } else {
        """
        others=${array}
        for i in \$(seq 1 \${#others[@]}); do
            somaticseq_split_vcf -infile \${others[i-1]} \\
                -snv "snvs_\${i}.vcf" -indel "indels_\${i}.vcf"
        done

        somaticseq_parallel.py \\
            paired \\
            ${args} \\
            --classifier-snv ${snp_classifier} \\
            --classifier-indel ${indel_classifier} \\
            --genome-reference ${reference} \\
            ${blacklist_flag} \\
            ${target_flag} \\
            --threads ${task.cpus} \\
            --tumor-bam-file ${tumor} \\
            --normal-bam-file ${normal} \\
            ${caller_flags}

        cp .command.out .log
        """
    }
    //
}
