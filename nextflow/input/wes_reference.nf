refdir = "/data/project/stemcell/shannc/reference"
params.ref = [
    "genome":"$refdir/genomes/GCF_000001405.40_GRCh38.p14_genomic.fna",
    "homopolymers_microsatellites": "$refdir/tool_specific/msihomopoly-GCF_000001405.40_GRCh38.p14_genomic.tsv",
    "delly_exclude": "TODO",
    "known_variants": ["$refdir/variants/dbSNP-GCF_000001405.40.gz",
                       "$refdir/variants/ExAC.r1.sites.vep.vcf.gz"
    ],
    "pileup": "$refdir/variants/ExAC.r1.sites.vep.vcf.gz",
    "genome_blacklist": "$refdir/blacklists/ENCFF356LFX_renamed.bed",
    "panel_of_normals": "TODO",
    "snpEff_db": "GRCh38.p14",
    "known_variants_somatic": ["$refdir/"],
    "baits": "$refdir/exome_kits/SureSelectHumanAllExonV6Hg38/Covered.bed.gz",
    "baits_il": "$refdir/exome_kits/SureSelectHumanAllExonV6Hg38/Covered.interval_list",
    "targets_il": "$refdir/exome_kits/SureSelectHumanAllExonV6Hg38/Regions.interval_list",
    "targets": "$refdir/exome_kits/SureSelectHumanAllExonV6Hg38/Regions.bed.gz",
    "delly_mappability": "$refdir/tool_specific/delly_mappability_GRCh38.fna.gz"
]

params.small_qc = ["accepted_filters": ["PASS"],
                   "min_tumor_depth": 10,
                   "max_normal_depth": 10,
                   "min_VAF": 0.10]

// Can't calculate AF, AD, DP or any of the standard INFO/FORMAT headers for SVs
// so will filter only with PASS
params.sv_qc = ["accepted_filters": ["PASS"]]
