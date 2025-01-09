refdir = "/data/project/stemcell/shannc/reference"
kitdir = "$refdir/exome_kits/SureSelectHumanAllExonV6Hg38"

genome = "$refdir/genomes/GRCh38.p14_filtered.fasta"
genome_gff = "$refdir/genomes/GRCh38.p14_filtered.gff"
genome_sdf = "$refdir/genomes/GRCh38.p14_filtered.sdf"
civic_cache = "/data/home/shannc/.cache/civic.json"
pandrugs2_cache = "/data/home/shannc/.cache/pandrugs2.json"

homopolymers_microsatellites= "$refdir/tool_specific/mishomopoly_GRCh38_filtered.tsv"
cosmic_reference = "${projectDir}/config/cosmic_signatures_v3.4-2024-12-26.csv"

targets= "$kitdir/Regions.bed.gz"
baits= "$kitdir/Covered.bed.gz"
baits_il= "$kitdir/Covered.interval_list"
targets_il= "$kitdir/Regions.interval_list"
baits_unzipped= "$kitdir/Unzipped_covered.bed"
genome_blacklist= "$refdir/blacklists/ENCFF356LFX_renamed.bed"
dbsnp = "$refdir/variants/dbSNP_renamed_germline.vcf.gz"

gnomad_all = "$refdir/variants/gnomADv4.1.0_all.vcf.gz"
gnomad_all_biallelic = "$refdir/variants/gnomADv4.1.0_all.vcf.gz"
gnomad_subset ="$refdir/variants/gnomADv4.1.0_Exomes/random/gnomADv4.1_subset.vcf.gz"
gnomad_subset_biallelic = "$refdir/variants/gnomADv4.1.0_Exomes/biallelic/all.vcf.gz"

dbsnp= "$refdir/variants/dbSNP_renamed_germline.vcf.gz"
pileup= "$refdir/variants/gnomADv4.1.0_all_biallelic.vcf.gz"
cnv_reference= "$refdir/variants/structural/aggregated_cnv.tsv"
msi_reference= "$refdir/variants/structural/aggregated_msi.tsv"

clingen_dosage = "$refdir/therapeutics/Clingen-Dosage-Sensitivity-2024-12-16.csv"
clingen_gene = "$refdir/therapeutics/Clingen-Gene-Disease-Summary-2024-12-19.csv"
