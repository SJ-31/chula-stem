miniforge3 = "/data/project/stemcell/shannc"

process {
    withName: "VEP" {
        container = "/data/home/shannc/tools/vep.sif" // Was installed with singularity
        ext.species = "homo_sapiens"
        ext.cache = "/data/home/shannc/.cache/vep"
        ext.args = ["--clin_sig_allele 0", // Reports alleles with known clinical significance
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

    withName: "GRIDSS" {
        label = "mid_mem"
        conda = "${miniforge3}/envs/gridss"
        cpus = 8 // Default is 8
    }

    withName: "CNVKIT|CNVKIT_PREP" {
        conda = "${miniforge3}/envs/cnvkit"
    }

    withName: "MUSE2" {
        container = "/data/home/shannc/tools/MuSE.sif"
        ext.cores = 8
        cpus = 8
    }

    withName: "OCTOPUS" {
        container = "/data/home/shannc/tools/octopus.sif"
        ext.args = ["--variant-discovery-mode ILLUMINA",
                    "--max-genotypes 400", // Should improve runtime
                    "--min-expected-somatic-frequency 0.05",
                    "--min-credible-somatic-frequency 0.02"] // Try to call somatic mutations
                    // at 2% VAF
        cpus = 4
    }

    withName: "MULTIQC" {
        container = "/data/home/shannc/tools/multiqc.sif"
    }

    withName: "FACETS_PILEUP|FACETS" {
        container = "/data/home/shannc/tools/facets.sif"
    }

    withName: "FACETS" {
        ext.args = ["--purity-cval 1000", "--cval 500", "--everything"]
        ext.genome = "hg38"
    }

    withName: "CLAIRS_TO" {
        container = "/data/home/shannc/tools/clairs-to.sif"
        ext.platform = "ilmn_ssrs"
        cpus = 4
    }

    withName: "CLAIRS" {
        container = "/data/home/shannc/tools/clairs.sif"
        // ext.platform = "ilmn_ss"
        ext.platform = "ilmn_ssrs"
        cpus = 4
    }

    withName: "SIGPROFILERASSIGNMENT" {
        conda = "${miniforge3}/envs/SigProfiler"
        ext.volume = "/data/project/stemcell/shannc/reference/tool_specific/sigprofiler"
    }

    withName: "DEEPSOMATIC" {
        container = "/data/home/shannc/tools/deepsomatic.sif"
        ext.model = "WES"
    }

    withName: "CLASSIFY_CNV" {
        ext.genome = "hg38"
    }

    withName: "STRELKA2|MANTA" {
        conda = "${miniforge3}/envs/strelka_manta"
        ext.args = ["--exome"]
    }

    withName: "CONCAT_VCF" {
        ext.publish = false
    }

}
