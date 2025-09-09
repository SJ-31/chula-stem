// * Required params
// Routine for processing the *recal files produced previously
workDir = "/data/project/stemcell/shannc/output/nf-work/PDAC"
includeConfig "$projectDir/input/wes_tumor_only_tools.config.nf"
params.outdir = "/data/project/stemcell/shannc/output/PDAC"
if (params.debug) {
    params.input = "$projectDir/input/manifests/pdac_DEBUG.csv"
} else {
    params.input = "$projectDir/input/manifests/pdac_prev.csv"
}
params.logdir = "/data/project/stemcell/shannc/output/PDAC/log"
params.routine = "wes"
params.cohort = "prev_cohort"
params.tumor_only = true
params.source_name = "TOOL_SOURCE"

// Map for pipeline resources
// See routine-specific documentation to see which need to be included

includeConfig "${projectDir}/input/references.config.nf"
params.report_data = ["hospital": "Dr.Nipan CU",
                      "physician": "Assoc. Prof. Nipan Israsena, M.D., Ph.D.",
                      "diagnosis": "Pancreatic ductal adenocarcinoma",
                      "sample_type": "organoid",
                      "collection_date": null,
]
params.report_text = ["front_page_details": "${projectDir}/reports/wes_front_page_details.txt",
                      "title": "Whole-exome sequencing results",
                      "disclaimer": "${projectDir}/reports/wes_disclaimer.txt",
                      "front_page_note": "${projectDir}/reports/wes_front_page_note.txt"]

genome_chr = "${refdir}/genomes/prev_assembly/Homo_sapiens_assembly38.fasta"
targets_chr = "${kitdir}/S07604514_Regions_1.bed.gz"
baits_chr = "${kitdir}/S07604514_Covered_1.bed.gz"
baits_chr_unzipped = "${kitdir}/S07604514_Covered.bed"
genome_gff_chr = "${refdir}/prev_assembly/hg38.ensGene.gtf"

targets_il_chr = "${kitdir}/Regions_chr.interval_list"
baits_il_chr = "${kitdir}/Covered_chr.interval_list"

dbsnp_chr = "${refdir}/variants/dbSNP_renamed_germline_chr.vcf.gz"
gnomad_subset_chr = "${refdir}/variants/gnomADv4.1.0_Exomes/random/gnomADv4.1_subset_chr.vcf.gz"
gnomad_subset_biallelic_chr = "${refdir}/variants/gnomADv4.1.0_Exomes/biallelic/all_chr.vcf.gz"
homopolymers_microsatellites_chr = "${refdir}/tool_specific/msihomopoly-Homo_sapiens_assembly38.tsv"

params.ref = ["genome": genome_chr,
              "homopolymers_microsatellites": homopolymers_microsatellites_chr,
              "targets": targets_chr,
              "baits": baits_chr,
              "baits_il": baits_il_chr,
              "targets_il": targets_il_chr,
              "baits_unzipped": baits_chr_unzipped,
              "genome_blacklist": genome_blacklist,
              "known_variants": [dbsnp_chr, gnomad_subset_chr],
              "genome_gff": genome_gff_chr,
              "germline": gnomad_subset_biallelic_chr,
              "dbsnp": dbsnp_chr,
              "pileup": gnomad_subset_biallelic_chr,
              "cnv_reference": cnv_reference,
              "msi_reference": msi_reference,
              "clingen_dosage": clingen_dosage,
              "civic_cache": null,
              "pandrugs2_cache": null,
              "cosmic_reference": cosmic_reference,
              "clingen_gene": clingen_gene]

// * Variant calling parameters

// ** Quality control filters
params.small_qc = ["accepted_filters": ["PASS"],
                   "min_tumor_depth": 10,
                   "max_normal_depth": 10,
                   "canonical": true,
                   "min_vaf": 0.10]
params.sv_qc = ["accepted_filters": ["PASS"]]

// * TEMPORARY
// <2025-01-08 Wed> Having difficulties with the therapy db caching, so for now
// will produce report without that data
process {
    withName: "REPORT" {
        ext.misc = ["plot": false, "show_therapies": true]
        ext.extra_paths = ["therapy_data_cache": "/data/home/shannc/.cache/2025-01-13_therapy_data.parquet"]
    }

    withName: "VEP" {
        container = "/data/home/shannc/tools/vep.sif" // Was installed with singularity
        ext.species = "homo_sapiens"
        ext.cache = "/data/project/stemcell/shannc/.cache/vep"
        ext.args = ["--clin_sig_allele 0", // Reports alleles with known clinical significance
                    "--synonyms /data/project/stemcell/shannc/reference/rename/vep_chr_synonyms.tsv", // NOTE: [2025-06-18 Wed] required for compatibility with previous PDAC variant calls
                    // in CLIN_SIG field
                    "--check_existing",// Check for existence of known variants co-located
                    "--pubmed", // Report pubmed ids for known variants
                    "--overlaps", // Report proportion and length of transcript overlapped
                    // by SV
                    "--exclude_predicted",
                    "--gene_phenotype", // Indicate if gene is associated with phenotype
                    "--regulatory", // Look for overlaps with regulatory regions
                    "--canonical", // Indicates if transcript is canonical for the gene
                    "--offline"
                    // "--plugin StructuralVariantOverlap,file=???" // Retrieve information from
                    // overlapping user-specified structural variants
                    ]
    }
}
