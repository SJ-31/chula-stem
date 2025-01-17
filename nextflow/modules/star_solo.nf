process STAR_SOLO {
    ext version: "2.7.11b"
    label "big_mem"

    publishDir "${meta.out}", mode:"copy", saveAs: params.saveFn
    publishDir "${meta.log}", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), val(reads) // reads
    val(reference) // Path to star index diretory
    val(strandedness) // forward|reverse|unstranded
    val(count) // Boolean value for whether or not to count reads
    val(gtf) // Genome gtf or gff file for counting reads
    val(barcodes) // 10x Genomics barcode whitelist
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

    check1 = file("${meta.out}/${output}")
    check2 = file("${meta.out}/${chimeric_junction}")
    check3 = file("${meta.out}/${counts}")

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

    strandedness_flag = params.strandedness ? "--soloStrand ${strandedness.toLowerCase()}" : " "
    solo_flag_list = [
        "--soloType CB_UMI_Simple", // Is type "Droplet"
        "--soloCBwhitelist ${barcodes}",
        strandedness_flag,
    ]
    solo_flag = Utl.overrideArgs(solo_flag_list, task.ext.args)

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

    args = task.ext.args.join(" ")
    if (check1.exists()) {
        """
        ln -sr ${check1} .
        if [[ -e ${check2} ]]; then ln -sr ${check2} .; fi
        if [[ -e ${check3} ]]; then ln -sr ${check3} .; fi
        ln -sr ${meta.log}/star_solo.log .
        """
    } else {
        """
        STAR ${args } --genomeDir ${reference} \\
            --readFilesIn ${reads[0]} ${reads[1]} \\
            --runThreadN ${task.cpus} \\
            --outSAMtype BAM SortedByCoordinate \\
            ${solo_flag} \\
            ${count_flags} \\
            ${compression_flag} \\
            ${strandedness_flag} \\
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

        get_nextflow_log.bash star_solo.log
        """
    }
    //
}
