miniforge3 = "/data/project/stemcell/shannc"

process {
    withName: "STAR" {
        // Will override any of the args provided in the module, namely those
        // related to fusion detection
        ext.args = []
        cpus = 4
    }

    withName: "STAR_FUSION" {
        container = "/data/home/shannc/tools/star-fusion.sif"
    }

    withName: "HTSEQ_COUNT" {
        conda = "${miniforge3}/envs/htseq"
    }

    withName: "TXIMPORT" {
        ext.args = ["--ignore_transcript_version"]
    }
}
