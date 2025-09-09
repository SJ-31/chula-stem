library(tidyverse)
library(glue)
library(httr2)

# [2025-09-09 Tue]
# Script to extract alterations to protein sequence features from VCF files,
#   scoring them by consequence
# Currently compatible with annotated VCFs produced by Ensembl vep
#
# Expected columns in the VCF
# - Consequence: Variant consequences defined by VEP
# - SYMBOL: HGNC gene symbol
# - Protein_position: amino acid position of the variant
# - Feature: Ensembl gene transcript ID

BASE_URL <- "https://rest.uniprot.org/uniprotkb"

# Cached results
SEEN_IDS <- c()
FEATURE_MAP <- tibble()

## * Dictionaries

SEQUENCE_FEATURES <- list(
  molecule_processing = list(
    ft_init_met = "Initiator methionine", # Cleavage of the initiator methionine
    ft_signal = "Signal", # Sequence targeting proteins to the secretory pathway or periplasmic space
    ft_transit = "Transit peptide", # Extent of a transit peptide for organelle targeting
    ft_propep = "Propeptide", # Part of a protein that is cleaved during maturation or activation
    ft_chain = "Chain", # Extent of a polypeptide chain in the mature protein
    ft_peptide = "Peptide" # Extent of an active peptide in the mature protein
  ),
  regions = list(
    ft_topo_dom = "Topological domain", # Location of non-membrane regions of membrane-spanning proteins
    ft_transmem = "Transmembrane", # Extent of a membrane-spanning region
    ft_intramem = "Intramembrane", # Extent of a region located in a membrane without crossing it
    ft_domain = "Domain", # Position and type of each modular protein domain
    ft_repeat = "Repeat", # Positions of repeated sequence motifs or repeated domains
    ft_zn_fing = "Zinc finger", # Position(s) and type(s) of zinc fingers within the protein
    ft_dna_bind = "DNA binding", # Position and type of a DNA-binding domain
    ft_region = "Region", # Region of interest in the sequence
    ft_coiled = "Coiled coil", # Positions of regions of coiled coil within the protein
    ft_motif = "Motif", # Short (up to 20 amino acids) sequence motif of biological interest
    ft_compbias = "Compositional bias" # Region of compositional bias in the protein
  ),
  site = list(
    ft_act_site = "Active site", # Amino acid(s) directly involved in the activity of an enzyme
    ft_binding = "Binding site", # Binding site for any chemical group (co-enzyme, prosthetic group, etc.)
    ft_site = "Site" # Any interesting single amino acid site on the sequence
  ),
  amino_acid_modifications = list(
    ft_non_std = "Non-standard residue", # Occurrence of non-standard amino acids (selenocysteine and pyrrolysine) in the protein sequence
    ft_mod_res = "Modified residue", # Modified residues excluding lipids, glycans and protein cross-links
    ft_lipid = "Lipidation", # Covalently attached lipid group(s)
    ft_carbohyd = "Glycosylation", # Covalently attached glycan group(s)
    ft_disulfid = "Disulfide bond", # Cysteine residues participating in disulfide bonds
    ft_crosslnk = "Cross-link" # Residues participating in covalent linkage(s) between proteins
  ),
  natural_variations = list(
    ft_var_seq = "Alternative sequence", # Amino acid change(s) producing alternate protein isoforms
    ft_variant = "Natural variant" # Description of a natural variant of the protein
  ),
  experimental_info = list(
    ft_mutagen = "Mutagenesis", # Site which has been experimentally altered by mutagenesis
    ft_unsure = "Sequence uncertainty", # Regions of uncertainty in the sequence
    ft_conflict = "Sequence conflict", # Description of sequence discrepancies of unknown origin
    ft_non_cons = "Non-adjacent residues", # Indicates that two residues in a sequence are not consecutive
    ft_non_ter = "Non-terminal residue" # Indicates that a residue is not the terminal residue of the complete protein
  ),
  secondary_structure = list(
    ft_helix = "Helix", # Helical regions within the experimentally determined protein structure
    ft_turn = "Turn", # Turns within the experimentally determined protein structure
    ft_strand = "Beta strand" # Beta strand regions within the experimentally determined protein structure
  )
)

CONSEQUENCE_SCORES <- list(
  transcript_ablation = 3,
  splice_acceptor_variant = 2,
  splice_donor_variant = 2,
  stop_gained = 3,
  frameshift_variant = 3,
  stop_lost = 3,
  start_lost = 3,
  transcript_amplification = 3,
  feature_elongation = 3,
  feature_truncation = 3,
  inframe_insertion = 2,
  inframe_deletion = 2,
  missense_variant = 2,
  mature_miRNA_variant = 2,
  protein_altering_variant = 2,
  splice_donor_5th_base_variant = 2,
  splice_region_variant = 2,
  splice_donor_variant = 2,
  splice_donor_region_variant = 2,
  splice_polypyrimidine_tract_variant = 2,
  incomplete_terminal_codon_variant = 3,
  start_retained_variant = 1,
  stop_retained_variant = 1,
  synonymous_variant = 1,
  coding_sequence_variant = 2,
  "5_prime_UTR_variant" = 1,
  "3_prime_UTR_variant" = 1,
  non_coding_transcript_exon_variant = 1,
  intron_variant = 1,
  NMD_transcript_variant = 1,
  non_coding_transcript_variant = 1,
  coding_transcript_variant = 2,
  upstream_gene_variant = 1,
  downstream_gene_variant = 1,
  TFBS_ablation = 2,
  TFBS_amplification = 2,
  TF_binding_site_variant = 2,
  regulatory_region_ablation = 2,
  regulatory_region_amplification = 2,
  regulatory_region_variant = 2,
  intergenic_variant = 1,
  sequence_variant = 1
)


## * Helper functions

#' Look up sequence features for the given uniprot accession
#'
#' @description
#' @param features Only lookup these features. See https://www.uniprot.org/help/sequence_annotation/
get_features <- function(uniprot_acc, features = NULL) {
  if (is.null(features)) {
    features <- names(flatten(SEQUENCE_FEATURES))
  }
  url <- glue("{BASE_URL}/{uniprot_acc}")
  req <- request(url) |>
    req_headers(accept = "application/json") |>
    req_url_query(fields = paste0(features, collapse = ",")) |>
    req_throttle(capacity = 60, fill_time_s = 60) |>
    req_retry(max_tries = 3)
  response <- req |>
    req_error(is_error = \(x) FALSE) |>
    req_perform()
  resp_body_json(response)
}

get_uniprot_ids <- function(
    tb,
    map_file,
    uniprot_colname = "uniprotswissprot",
    ensembl_id_colname = "ensembl_transcript_id") {
  tmp <- read_tsv(map_file) |>
    rename(Feature = ensembl_id_colname) |>
    distinct()
  tmp <- tmp[!is.na(tmp[[uniprot_colname]]), ]
  tmp |> filter(Feature %in% tb$Feature)
}


make_feature_tb <- function(symbol, feature_dct) {
  feature_df <- lapply(feature_dct$features, \(x) {
    start <- x$location$start$value
    type <- str_to_lower(x$type) |> str_replace_all(" ", "_")
    tibble(
      symbol = symbol,
      name = glue("{symbol}.{type}.{start}"),
      start = start,
      end = x$location$end$value
    )
  }) |>
    bind_rows()
}

#' Retrieve multiple lists of sequence features
#'
#' @description
#' @param feature_tb gene tb (unique) with at least two columns:
#'  1. uniprot accessions and 2. symbol names
#'
get_all_features <- function(
    feature_tb,
    features = NULL,
    uniprot_colname = "uniprotswissprot",
    symbol_colname = "hgnc_symbol") {
  tb <- apply(feature_tb, 1, \(row) {
    if (!row[uniprot_colname] %in% SEEN_IDS) {
      SEEN_IDS <<- c(row[uniprot_colname], SEEN_IDS)
      lookup <- get_features(row[uniprot_colname], features)
      if (!"features" %in% lookup) {
        make_feature_tb(row[symbol_colname], lookup)
      } else {
        NULL
      }
    } else {
      NULL
    }
  }) |>
    bind_rows()
  tb
}


add_feature_score <- function(tb, how = "consequence") {
  if (how == "consequence") {
    mutate(
      tb,
      feature_score = map_dbl(Consequence, \(csq) {
        if (is.na(csq)) {
          return(0)
        }
        if (str_detect(csq, ";")) {
          consequences <- str_split_1(csq, ";")
          map_dbl(consequences, \(c) CONSEQUENCE_SCORES[[c]]) |>
            max()
        } else {
          lookup <- CONSEQUENCE_SCORES[[csq]]
          if (is.null(lookup)) {
            0
          } else {
            lookup
          }
        }
      })
    )
  } else {
    mutate(tb, feature_score = 1)
  }
}

## * Main routine

#' Encode gene sequence attributes in `sample_tb` as features, whose values are determined
#'  by their mutations in the sample
#'
#' @description
#' @param sample_tb annotated VCF tibble from Ensembl vep
#' @param only_symbols Only encode the given symbols. By default, tries to encode all
#' @param how encoding method, see `add_feature_score`
encode_w_uniprot <- function(
    sample_tb,
    feature_map,
    how = "consequence") {
  extra <- list()
  joined <- lapply(unique(sample_tb$SYMBOL), \(sym) {
    sample_filtered <- sample_tb |> filter(SYMBOL == sym)
    ftb_filtered <- feature_map |> filter(symbol == sym)
    intersected <- left_join(
      sample_filtered,
      ftb_filtered,
      by = join_by(
        between(x$Protein_position, y_lower = y$start, y_upper = y$end)
      )
    )
    noncoding <- sample_filtered |>
      filter(is.na(Protein_position)) |>
      add_feature_score(how)
    other <- sample_filtered |>
      filter(!Loc %in% intersected$Loc & !Loc %in% noncoding$Loc) |>
      add_feature_score(how)
    extra[[glue("{sym}.other_protein_coding")]] <<- sum(
      other$feature_score
    )
    extra[[glue("{sym}.noncoding")]] <<- sum(noncoding$feature_score)
    intersected
  }) |>
    bind_rows() |>
    filter(!is.na(name))
  encoded <- joined |>
    select(name, Consequence) |>
    add_feature_score() |>
    select(-Consequence) |>
    group_by(name) |>
    summarise(feature_score = sum(feature_score)) |>
    pivot_wider(names_from = name, values_from = feature_score)
  bind_cols(as_tibble(extra), encoded)
}

encode_w_vep <- function(sample_tb, allowed_consequence = NULL) {
  if (!is.null(allowed_consequence)) {
    mask <- map_lgl(sample_tb$Consequence, \(c) {
      if (str_detect(c, ";")) {
        splits <- str_split_1(c, ";")
        length(intersect(splits, allowed_consequence)) >= 1
      } else {
        c %in% allowed_consequence
      }
    })
    sample_tb <- sample_tb[mask, ]
  }
}

encode_multiple_samples <- function(
    method,
    sample_path_map,
    map_file,
    features = NULL,
    only_symbols = NULL,
    allowed_consequence = NULL,
    uniprot_colname = "uniprotswissprot",
    ensembl_id_colname = "ensembl_transcript_id") {
  lapply(names(sample_path_map), \(s) {
    sample_tb <- read_tsv(sample_path_map[[s]]) |>
      mutate(Protein_position = as.numeric(Protein_position))
    if (!is.null(only_symbols)) {
      sample_tb <- filter(sample_tb, SYMBOL %in% only_symbols)
    }
    if (method == "uniprot_feature") {
      cur_fmap <- get_uniprot_ids(
        filter(sample_tb, !is.na(Protein_position) & !is.na(SYMBOL)),
        map_file = map_file,
        uniprot_colname = uniprot_colname,
        ensembl_id_colname = ensembl_id_colname
      )
      FEATURE_MAP <<- bind_rows(
        list(FEATURE_MAP, get_all_features(cur_fmap, features = features))
      )
      encode_w_uniprot(sample_tb, feature_map = FEATURE_MAP) |>
        mutate(sample = s)
    } else if (method == "consequence") {
      encode_w_vep(sample_tb, allowed_consequence = allowed_consequence)
    } else {
      stop("Given encoding method not supported!")
    }
  }) |>
    bind_rows() |>
    mutate(across(where(is.numeric), \(x) replace_na(x, 0))) |>
    relocate(sample, .before = everything())
}

## * Entry point

from_config <- function(config_file, output) {
  if (str_detect(config_file, "\\.yaml$")) {
    config <- yaml::read_yaml(config_file)
  } else {
    config <- jsonlite::read_json(config_file)
  }
  if (is.null(config$features)) {
    features <- c(
      names(discard_at(SEQUENCE_FEATURES$regions, \(x) x == "ft_compbias")),
      names(SEQUENCE_FEATURES$site),
      names(SEQUENCE_FEATURES$secondary_structure),
      "ft_init_met",
      "ft_signal",
      "ft_transit",
      "ft_chain"
    )
  } else {
    features <- config$features
  }
  result <- encode_multiple_samples(
    method = config$encoding_method,
    sample_path_map = config$samples,
    map_file = config$transcript_id2uniprot$file,
    uniprot_colname = config$transcript_id2uniprot$uniprot_colname,
    ensembl_id_colname = config$transcript_id2uniprot$ensembl_id_colname,
    allowed_consequence = config$allowed_consequence,
    only_symbols = config$symbols,
    features = features
  )
  write_tsv(result, output)
}

if (sys.nframe() == 0) {
  library(optparse)
  parser <- OptionParser()
  parser <- add_option(
    parser,
    c("-i", "--input"),
    type = "character",
    help = "Config file input"
  )
  parser <- add_option(
    parser,
    c("-e", "--encoding_method"),
    type = "character",
    help = "How to encode variants into features. One of 'uniprot_feature', 'vep_consequence'",
    default = "uniprot_feature"
  )
  parser <- add_option(
    parser,
    c("-o", "--output"),
    type = "character",
    help = "Output file name"
  )
  args <- parse_args(parser)
  from_config(args$input, args$output)
}
