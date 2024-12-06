t = "/data/project/stemcell/shannc/tests"
refdir = "/data/project/stemcell/shannc/reference"
params.tools = "$projectDir/input/wes_tools.config.nf"
params.outdir = "$t/wes_test2"
params.logdir = "$t/wes_test2/log"
params.ref = ["genome": "$t/wes_test2/GRCh38.p14_filtered.fna",
              "known_variants": ["$refdir/variants/dbSNP_renamed.vcf.gz"],
              "delly_exclude": "$refdir/tool_specific/hg38.excl.tsv",
              "mappability": "$refdir/tool_specific/delly_mappability_GRCh38.fna.gz",
              "homopolymers_microsatellites": "$refdir/tool_specific/mishomopoly_GRCh38_filtered.tsv",
              "germline": "/data/project/stemcell/shannc/tests/wes_test2/test_variant/test.vcf.gz",
              "pileup": "/data/project/stemcell/shannc/tests/wes_test2/test_variant/test.vcf.gz", // Only chr1 for speed
              "dbsnp": "$refdir/variants/dbSNP_renamed.vcf.gz",
              "genome_blacklist": "$refdir/blacklists/ENCFF356LFX_renamed.bed",
              "targets": "$refdir/exome_kits/SureSelectHumanAllExonV6Hg38/Regions.bed.gz",
              "baits": "$refdir/exome_kits/SureSelectHumanAllExonV6Hg38/Covered.bed.gz",
            "baits_il": "$refdir/exome_kits/SureSelectHumanAllExonV6Hg38/Covered.interval_list",
            "targets_il": "$refdir/exome_kits/SureSelectHumanAllExonV6Hg38/Regions.interval_list",
              "baits_unzipped": "$refdir/exome_kits/SureSelectHumanAllExonV6Hg38/Unzipped_covered.bed"
]
params.routine = "wes"
params.input = "/data/project/stemcell/shannc/tests/wes_test2/manifest.csv"

params.small_qc = ["accepted_filters": ["PASS"],
                   "min_tumor_depth": 10,
                   "max_normal_depth": 10,
                   "min_VAF": 0.10]

// Can't calculate AF, AD, DP or any of the standard INFO/FORMAT headers for SVs
// so will filter only with PASS
params.sv_qc = ["accepted_filters": ["PASS"]]
