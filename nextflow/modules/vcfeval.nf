process VCFEVAL {
    ext version: "3.12.1"

    publishDir "${meta.out}", mode:"copy", saveAs: params.saveFn
    publishDir "${meta.log}", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), path(baseline), path(query)
    val(reference) // Genome reference in SDF format
    val(regions)
    val(module_number)
    //

    output:
    path("*.log")
    //

    script:
    output = Utl.getName(module_number, meta, "vcfeval")
    check = file("${meta.out}/${output}")
    args = task.ext.args.join(" ")
    regions_flag = regions ? " --bed-regions=${regions} " : ""
    if (task.ext.rename_chrs) {
    // For when baseline and query use different chromosome naming conventions
        rename_file = task.ext.rename_chrs.file
          // File of "old_name new_name\n" pairs
        target = task.ext.rename_chrs.rename_target == "base" ? baseline : query
        chr_rename = "bcftools annotate --rename-chrs ${rename_file} ${target} -O z > temp1.vcf.gz"
        if (target == query) {
            query = "temp1.vcf.gz"
        } else {
            baseline = "temp1.vcf.gz"
        }
    } else {
        chr_rename = ""
    }

    if (task.ext.rename_samples) {
        // For when the baseline and query use different sample naming conventions
        rename_file = task.ext.rename_samples.file
        if (task.ext.rename_samples.rename_target == "both") {
            rename_query = "bcftools reheader --samples ${rename_file} ${query} -O z > query.vcf.gz"
            rename_base = "bcftools reheader --samples ${rename_file} ${baseline} -O z > base.vcf.gz"
            sample_rename = "${rename_query}; ${rename_base}"
            query = "query.vcf.gz"
            baseline = "base.vcf.gz"
        } else {
            target = task.ext.rename_samples.rename_target == "base" ? baseline : query
            sample_rename = "bcftools reheader --samples ${rename_file} ${target} -O z > temp2.vcf.gz"
            if (target == query) {
                query = "temp2.vcf.gz"
            } else {
                baseline = "temp2.vcf.gz"
            }
        }
    } else {
        sample_rename = ""
    }
    if (check.exists()) {
        """
        cp -r ${check} .
        ln -sr ${meta.log}/vcfeval.log .
        """
    } else {
        """
        ${chr_rename}

        ${sample_rename}

        rtg vcfeval \\
            ${args} \\
            ${regions_flag} \\
            --output=${output} \\
            --baseline=${baseline} \\
            --calls=${query} \\
            --template=${reference}

        get_nextflow_log.bash vcfeval.log
        """
    }
    //
}
