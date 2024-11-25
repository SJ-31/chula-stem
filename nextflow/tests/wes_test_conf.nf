t = "/data/project/stemcell/shannc/tests"
refdir = "/data/project/stemcell/shannc/reference"
params.tools = "$projectDir/input/wes_tools.nf"
params.outdir = "$t/wes_test2"
params.logdir = "$t/wes_test2/log"
params.ref = ["genome": "$t/wes_test2/GRCh38.p14_filtered.fna",
              "known_variants": ["$refdir/variants/dbSNP_renamed.vcf.gz"],
              "delly_exclude": "$refdir/tool_specific/hg38.excl.tsv",
              "mappability": "$refdir/tool_specific/delly_mappability_GRCh38.fna.gz",
              "homopolymers_microsatellites": "$refdir/tool_specific/mishomopoly_GRCh38_filtered.tsv",
              "pileup": "/data/project/stemcell/shannc/tests/wes_test2/test_variant/test.vcf.gz",
              "dbsnp": "$refdir/variants/dbSNP_renamed.vcf.gz",
              "genome_blacklist": "$refdir/blacklists/ENCFF356LFX_renamed.bed"
]
params.routine = "wes"
params.input = "/data/project/stemcell/shannc/tests/wes_test2/manifest.csv"
