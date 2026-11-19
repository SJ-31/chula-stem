include { OPTITYPE } from "./modules/optitype.nf"
include { MERGE_FASTQ } from "./modules/merge_fastq.nf"
include { PURECN_COMPLETE } from "./subworkflows/purecn_complete.nf"
include { SAREK_MANIFEST } from "./modules/sarek_manifest.nf"
include { CONCAT_VCF } from "./modules/concat_vcf.nf"
include { SIGPROFILERASSIGNMENT ; SIGPROFILERASSIGNMENT_COLLECT } from "./modules/sigprofilerassignment.nf"
include { EXTRACT_FORMAT_FIELDS } from "./modules/extract_format_fields.nf"
include { UPDATE_GENE_CACHE } from "./modules/update_gene_cache.nf"
include { MERGE_CLASSIFY_CNV } from "./modules/merge_classify_cnv.nf"
include { FORMAT_VEP_VCF } from "./modules/format_vep_vcf.nf"

nextflow.enable.types = true
params {
    sarek_config: Path
    defaults: Map
    update_cache: Boolean
    genome: Path
    sigprofiler_volume: Path
    sigprofiler_exclude: Path?
    gene_data_cache: String
    open_targets: String?
    panel_of_normals: Path?
    baits: Path?
    mappability: Path
    purecn_bait_intervals: Path?
    purecn_mapping_bias: Path?
    purecn_normaldb: Path?
    snp_blacklist: Path?
}

// NOTE: https://github.com/nf-core/sarek/blob/master/conf/igenomes.config for the igenomes references


workflow {

    main:
    sarek_out = SAREK_MANIFEST(params.sarek_config)

    preprocessed = sarek_out
        .flatMap { v ->
            v.preprocessed.splitCsv(header: true)
        }
        .map { row ->
            record(
                id: row.id,
                status: row.status as Integer,
                recal: row.recal,
                recal_index: row.recal_index,
                fastq_1: row.fastq_1,
                fastq_2: row.fastq_2,
                split_fastq: row.split_fastq,
            )
        }

    tool_outs = sarek_out.flatMap { v ->
        v.tool_outs.splitCsv(header: true)
    }

    merged = MERGE_FASTQ(
        preprocessed.filter { v -> v.split_fastq }.map { v ->
            record(
                id: v.id,
                fastq_dir: v.fastq_1,
                status: v.status,
            )
        }
    )

    fastq = preprocessed
        .filter { v -> !v.split_fastq }
        .map { v ->
            record(
                id: v.id,
                fastq_1: file(v.fastq_1),
                fastq_2: file(v.fastq_2),
                status: v.status,
            )
        }
        .mix(merged)

    optitype = OPTITYPE(fastq.filter { v -> v.status == 0 }, "dna")

    to_purecn = tool_outs
        .map { row ->
            record(id: row.id, mutect2_vcf: row.mutect2_unfiltered)
        }
        .join(
            preprocessed.filter { v -> v.status == 0 }.map { v ->
                record(
                    id: v.id,
                    normal_bam: v.recal,
                    normal_bam_index: v.recal_index,
                )
            },
            by: "id"
        )
        .join(
            preprocessed.filter { v -> v.status == 1 }.map { v ->
                record(
                    id: v.id,
                    tumor_bam: v.recal,
                    tumor_bam_index: v.recal_index,
                )
            },
            by: "id"
        )

    // TODO: check if the GEM mappability file used by Control-FREEC
    // also works with PureCN

    // purecn_complete = PURECN_COMPLETE(
    //     samples: to_purecn,
    //     panel_of_normals: params.ref.panel_of_normals,
    //     genome: params.ref.genome,
    //     baits: params.ref.baits,
    //     mappability: params.ref.mappability,
    //     purecn_bait_intervals: params.ref.purecn_bait_intervals,
    //     purecn_mapping_bias: params.ref.purecn_mapping_bias,
    //     purecn_normaldb: params.ref.purecn_normaldb,
    //     snp_blacklist: params.ref.snp_blacklist,
    //     default_purity: params.defaults.purity,
    //     default_ploidy: params.defaults.ploidy
    // )
    
    sv_tumor_to_concat = tool_outs
        .filter { v -> v.output_type != "germline" }
        .filter { v -> v.vep_manta_somatic || v.vep_tiddit }
        .map { v ->
            record(
                id: v.id,
                vcfs: ([v?.vep_manta_somatic, v?.vep_tiddit].findAll { el -> el }.toList()),
            )
        }

    sv_tumor_concat = CONCAT_VCF(sv_tumor_to_concat, "sv_")

    to_extract_format = tool_outs.map { v ->
        def res = [:]
        ["deepvariant", "freebayes", "strelka", "tiddit", "deepsomatic", "manta_somatic", "muse", "mutect2", "strelka_indels", "haplotypecaller", "deepvariant"].each { el ->
            if (v[el] && file(v[el]).exists()) {
                res[el] = file(v[el])
            }
        }
        record(
            id: v.id,
            caller2vcfs: res,
            suffix: "_${v.output_type}",
            normal_sample: v.normal_sample,
            tumor_sample: v.tumor_sample,
        )
    }
    format_fields = EXTRACT_FORMAT_FIELDS(to_extract_format)

    cnv_to_merge = tool_outs
        .map { v ->
            def res = [:]
            ["purecn_cns", "controlfreec_cns", "ascat_cns", "cnvkit_cns"].each { el ->
                if (v[el] && file(v[el]).exists()) {
                    res[el.replace("_cns", "")] = file(v[el])
                }
            }
            record(id: v.id, caller2cnv_file: res)
        }
        .filter { v -> v.caller2cnv_file.size() > 1 }
    // TODO: [2026-11-17 Tue] put this into the report
    merged_cnvs = MERGE_CLASSIFY_CNV(cnv_to_merge, "hg38")


    to_vep_tsv = tool_outs
        .filter { v -> v.vep_consensus }
        .map { v -> record(id: v.id, vcf: file(v.vep_consensus)) }
        .mix(
            sv_tumor_concat.map { el -> record(id: el.id, vcf: el.concat) }
        )
    vep_tsvs = FORMAT_VEP_VCF(
        to_vep_tsv,
        ["CALLERS", "NCALLERS", "SVTYPE"],
        "CSQ",
    )
    
    mutational_signatures = SIGPROFILERASSIGNMENT(
        tool_outs
            .filter { v -> v.output_type != "germline" }
            .map { v -> record(id: v.id,
                               vcfs: [file(v.consensus)]) },
        true,
        params.sigprofiler_exclude,
        params.sigprofiler_volume.toUriString(),
    )

    if (params.update_cache) {
        cache_proc = UPDATE_GENE_CACHE(vep_tsvs.map { v -> v.tsv }.collect(),
                                       params.ref.gene_data_cache,
                                       params.ref?.open_targets,
                                       "",
                                       true
        )
        report_cache = cache_proc.map { v -> v.cache } 
    } else {
        report_cache = channel.value(params.ref.gene_data_cache)
    }
    
    
    publish:


    // combined_signatures = SIGPROFILERASSIGNMENT_COLLECT(
    //     mutational_signatures.map { v -> v.activities }.collect(),
    //     "TODO", // Put the sample pattern here once you figure out what it is
    // )
    optitype = optitype
    format_fields = format_fields
    vep_tsvs = vep_tsvs
    sarek_out = sarek_out
    mutational_signatures = mutational_signatures
    merged_cnvs = merged_cnvs
}

output {
    mutational_signatures {
        path { v -> v.dir >> "sigprofilerassignment/${v.id}/${v.dir.name}" 
        }
    }
    merged_cnvs {
        path { v ->
            v.classified >> "merged_cnvs/${v.id}/${v.classified.name}"
            v.concordant >> "merged_cnvs/${v.id}/${v.concordant.name}"
            v.discordant >> "merged_cnvs/${v.id}/${v.discordant.name}"
        }
    }
    sarek_out {
        path { v ->
            v.tool_outs >> v.tool_outs.name
            v.preprocessed >> v.preprocessed.name
        }
    }
    optitype {
        path { v ->
            v.optitype_csv >> "optitype/${v.id}/${v.optitype_csv.name}"
            v.optitype_pdf >> "optitype/${v.id}/${v.optitype_pdf.name}"
        }
    }
    vep_tsvs {
        path { v -> v.tsv >> "annotation/${v.id}/${v.tsv.name}" }
    }
    format_fields {
        path { v -> v.tsv >> "format_fields/${v.id}/${v.tsv.name}" }
    }
}
