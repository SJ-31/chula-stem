process STAR {
    ext version: "2.7.11b"
    label "big_mem"

    publishDir "${meta.out}", mode:"copy", saveAs: params.saveFn
    publishDir "${meta.log}", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), val(reads)
    val(reference) // Path to star index diretory
    val(strandedness) // forward|reverse|unstranded
    val(mark_duplicates) // Mode for marking duplicates (see manual)
    // "-" for none, "UniqueIdentical", "UniqueIdenticalNotMulti"
    val(count) // Boolean value for whether or not to count reads
    val(gtf) // Genome gtf or gff file for counting reads
    val(module_number)
    //

    output:
    tuple val(meta), path(output), emit: mapped
    tuple val(meta), path(counts), emit: counts, optional: true
    tuple val(meta), path(chimeric_junction), emit: chimeric, optional: true
    path("*.log")
    //

    script:
    output = Utl.getName(module_number, meta, null, "bam")

    j_suffix = params.detect_fusion ? "out.junction" : "empty"
    c_suffix = count && gtf ? "tsv" : "empty"

    chimeric_junction = Utl.getName(module_number, meta, "STAR_Chimeric", j_suffix)
    counts = Utl.getName(module_number, meta, "STAR_Counts", c_suffix)

    check = file("${meta.out}/${output}")

    if (count && gtf) {
        count_flags = " --quantMode GeneCounts --sjdbGTFfile ${gtf}"
        if (strandedness == "forward") {
            get_counts_command = "cut -f 1,3 ReadsPerGene.out.tab > ${counts}"
        } else if (strandedness == "reverse") {
            get_counts_command = "cut -f 1,4 ReadsPerGene.out.tab > ${counts}"
        } else {
            get_counts_command = "cut -f 1,2 ReadsPerGene.out.tab > ${counts}"
        }
    } else {
        count_flags = ""
        get_counts_command = ""
    }

    fusion_flag_list = [
            "--chimSegmentMin 12",
            "--chimJunctionOverhangMin 8",
            "--chimOutJunctionFormat 1",
            "--alignSJDBoverhangMin 10",
            "--alignMatesGapMax 100000",
            "--alignIntronMax 100000",
            "--alignSJstitchMismatchNmax 5 -1 5 5",
            "--outSAMattrRGline ID:GRPundef",
            "--chimMultimapScoreRange 3",
            "--chimScoreJunctionNonGTAG -4",
            "--chimMultimapNmax 20",
            "--chimNonchimScoreDropMin 10",
            "--peOverlapNbasesMin 12",
            "--peOverlapMMp 0.1",
            "--alignInsertionFlush Right",
            "--alignSplicedMateMapLminOverLmate 0",
            "--alignSplicedMateMapLmin 30",
    ]
    fusion_flag = params.detect_fusion ? Utl.overrideArgs(fusion_flag_list, task.ext.args) : ""
    compression_flag = params.fastq_uncompress ? " --readFilesCommand ${params.fastq_uncompress}" : ""
    duplicates_flag = Utl.overrideArgs(["--bamRemoveDuplicatesType ${mark_duplicates}"],
                                       task.ext.args)
    args = task.ext.args.join(" ")
    if (check.exists()) {
        """
        ln -sr ${check} .
        ln -sr ${meta.log}/star.log .
        """
    } else {
        """
        STAR ${args } --genomeDir ${reference} \\
            --readFilesIn ${reads[0]} ${reads[1]} \\
            --runThreadN ${task.cpus} \\
            --outSAMtype BAM SortedByCoordinate \\
            ${count_flags} \\
            ${compression_flag} \\
            ${duplicates_flag} \\
            ${fusion_flag}

        if [[ ${params.detect_fusion} == "true" ]]; then
            mv Chimeric.out.junction ${chimeric_junction}
        fi

        ${get_counts_command}

        gatk AddOrReplaceReadGroups \\
            -I Aligned.sortedByCoord.out.bam \\
            -O ${output} \\
            --RGLB $meta.RGLB \\
            --RGPL $meta.RGPL \\
            --RGPU $meta.RGPU \\
            --RGSM $meta.RGSM

        get_nextflow_log.bash star.log
        """
    }
    //
}
