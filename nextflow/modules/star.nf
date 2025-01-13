process STAR {
    ext version: "2.7.11b"
    label "big_mem"

    publishDir "${meta.out}", mode:"copy", saveAs: params.saveFn
    publishDir "${meta.log}", mode: "copy", pattern: "*.log"

    input:
    tuple val(meta), val(reads)
    val(reference) // Path to star index diretory
    val(strandedness)
    val(module_number)
    //

    output:
    tuple val(meta), path(output), emit: mapped
    tuple val(meta), path(chimeric_junction), emit: chimeric, optional: true
    path("*.log")
    //

    script:
    output = Utl.getName(module_number, meta, null, "bam")

    j_suffix = params.detect_fusion ? "out.junction" : "empty"
    chimeric_junction = Utl.getName(module_number, meta, "STAR_CHIMERIC", j_suffix)
    strandedness_flag = params.strandedness ? "--soloStrand ${strandedness.toLowerCase()}" : " "

    check1 = file("${meta.out}/${output}")
    check2 = file("${meta.out}/${chimeric_junction}")

    other_args = task.ext.args.collect({ it.split()[0] })
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
    ].findAll({ !other_args.contains(it.split()[0]) })
    fusion_flag = params.detect_fusion ? fusion_flag_list.join(" ") : ""
    compression_flag = params.fastq_uncompress ? " --readFilesCommand ${params.fastq_uncompress}" : ""

    args = task.ext.args.join(" ")
    if (check1.exists() && check2.exists()) {
        """
        ln -sr ${check1} .
        ln -sr ${check2} .
        ln -sr ${meta.log}/star.log .
        """
    } else {
        """
        STAR ${args } --genomeDir ${reference} \\
            --readFilesIn ${reads[0]} ${reads[1]} \\
            --runThreadN ${task.cpus} \\
            --outSAMtype BAM SortedByCoordinate \\
            ${compression_flag} \\
            ${strandedness_flag} \\
            ${fusion_flag}

        if [[ ${params.detect_fusion} == "true" ]]; then
            mv Chimeric.out.junction ${chimeric_junction}
        else
            touch ${chimeric_junction}
        fi

        mv Aligned.sortedByCoord.out.bam ${output}

        get_nextflow_log.bash star.log
        """
    }
    //
}
