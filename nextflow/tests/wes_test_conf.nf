t = "/data/project/stemcell/shannc/tests"
params.tools = "$projectDir/input/wes_tools.nf"
params.outdir = "$t/wes_test2"
params.logdir = "$t/wes_test2/log"
params.ref = ["genome": "$t/wes_test2/GRCh38.p14_filtered.fna",
        "known_variants": ["/data/project/stemcell/shannc/reference/variants/dbSNP_renamed.vcf.gz"]
]
params.routine = "wes"
params.input = "/data/project/stemcell/shannc/tests/wes_test2/manifest.csv"
