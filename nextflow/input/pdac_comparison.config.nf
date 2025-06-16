// 2024-12-30 Config for comparing the small variant calls between PDAC runs

params.outdir = "/data/project/stemcell/shannc/output/PDAC_comparison"
params.input = ""
params.routine = "compare_variants"

includeConfig "${projectDir}/input/references.config.nf"

params.ref = ["genome_sdf": genome_sdf, "targets": targets]

process {
    withName: "VCFEVAL" {
        ext.args = ["--decompose",
                    "--output-mode=combine",
                    "--squash-ploidy" // Recommended for comparing somatic variants due
                    // to their variable ploidy levels
        ]
        ext.rename_samples = ["file": "/data/project/stemcell/shannc/reference/rename/remove_chr.tsv"
                              , rename_target: "base"]
        ext.rename_chrs = ["file": "${projectDir}/config/sample_rename.txt",
                           rename_target: "both"]
    }
}
