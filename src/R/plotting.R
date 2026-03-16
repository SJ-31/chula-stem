library(tidyverse)
library(checkmate)
library(patchwork)
library(pointblank)
library(glue)
library(grid)
library(ggplot2)
R_SRC <- Sys.getenv("R_SRC")
UTILS_FILE <- paste0(R_SRC, "/", "utils.R")
if (file.exists(UTILS_FILE)) {
  source(UTILS_FILE)
}
PY_SRC <- Sys.getenv("PY_SRC")


#' Returns "default" if the value of "key" in "list" is NULL
#'
lget <- function(list, key, default) {
  val <- list[[key]]
  if (is.null(val)) {
    default
  } else {
    val
  }
}

format_chr_data <- function(tb, divide_by = 1000) {
  chr_levels <- c(as.character(1:22), "X", "Y") %>% paste0("chr", .)
  if (!str_detect(tb$chromosome[1], "^chr")) {
    tb$chromosome <- paste0("chr", tb$chromosome)
  }
  tb |>
    mutate(
      start = start / divide_by,
      end = end / divide_by,
      mid = (end - start) / 2 + start,
      chromosome = factor(tb$chromosome, levels = chr_levels)
    )
}

plot_cnvkit_gene <- function(cnv_tb, locus) {}

parse_chr_range <- function(chr_range) {
  splits <- str_split_1(chr_range, ":")
  chr <- splits[1]
  if (length(splits) > 1) {
    endpoints <- str_split_1(splits[2], "-")
    list(
      chr = chr,
      ir = IRanges(as.numeric(endpoints[1]), as.numeric(endpoints[2]))
    )
  } else {
    list(chr = chr, ir = NULL)
  }
}

plot_cnvkit <- function(cns, cnr, chr, sizing, output) {
  cn_calls <- read_tsv(cns) |> format_chr_data()
  ratios <- read_tsv(cnr) |> format_chr_data()
  range_given <- FALSE
  if (!is.null(chr) && length(chr) >= 1) {
    if (str_detect(chr, "-")) range_given <- TRUE
    range <- parse_chr_range(chr)
    cn_calls <- cn_calls |> dplyr::filter(chromosome %in% range$chr)
    ratios <- ratios |> dplyr::filter(chromosome %in% range$chr)
  }
  merged <- ratios |>
    inner_join(
      cn_calls,
      by = join_by(chromosome, between(mid, y$start, y$end)),
      suffix = c("", "_right")
    )
  if (range_given) {
    vals <- map2(merged$start_right, merged$end_right, \(s, e) {
      intersect <- pintersect(
        IRanges(s, e),
        range$ir,
        resolve.empty = "start.x"
      )
      if (width(intersect) > 0) {
        tibble(start_right = start(intersect), end_right = end(intersect))
      } else {
        tibble(start_right = NA, end_right = NA)
      }
    }) |>
      bind_rows()
    merged <- merged |>
      dplyr::select(-start_right, -end_right) |>
      bind_cols(vals) |>
      dplyr::filter(start >= start(range$ir) & end <= end(range$ir))
  }
  plot <- merged |>
    ggplot(aes(x = mid, y = log2, alpha = weight), stat = "identity") +
    ylab("Copy Number Ratio (log2)") +
    geom_point(size = 0.5) +
    scale_y_continuous(
      breaks = scales::breaks_pretty(n = 15),
      minor_breaks = scales::minor_breaks_n(5),
      limits = c(min(ratios$log2), max(ratios$log2)),
    ) +
    geom_segment(
      aes(
        x = start_right,
        xend = end_right,
        y = log2_right,
        yend = log2_right,
        color = factor(cn)
      ),
      alpha = 0.5,
      linewidth = 1.5
    ) +
    guides(
      color = guide_legend(title = "Integer Copy Number"),
      alpha = guide_legend(title = "Weight")
    ) +
    xlab(NULL)
  if (is.null(chr) || length(chr) > 1) {
    without_chr <- str_replace(chr, "^chr", "")
    title <- ifelse(
      missing(chr),
      "All Chromosomes",
      glue("Chromosomes {paste(without_chr, collapse = ', ')}")
    )
    plot <- plot +
      facet_wrap(
        vars(chromosome),
        strip.position = "bottom",
        scales = "free_x",
        ncol = lget(sizing, "ncol", 3)
      ) +
      theme_bw() +
      theme(
        text = element_text(size = 15),
        axis.line = element_line(colour = "grey"),
        strip.placement = "outside",
        panel.grid.minor.x = element_blank(),
        legend.title = element_text(face = "bold"),
        title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        panel.grid.major.x = element_blank(),
        axis.ticks.y = element_line(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.background = element_rect(fill = "white", colour = "white")
      ) +
      labs(title = title)
  } else {
    without_chr <- str_replace(chr, "^chr", "")
    plot <- plot +
      xlab("Position (mb)") +
      labs(title = glue("Chromosome {without_chr}"))
  }
  if (range_given) {
    plot <- plot + scale_x_continuous(limits = c(range$start, range$end))
  }
  if (!is.null(output)) {
    ggsave(
      output,
      plot,
      dpi = lget(sizing, "dpi", 300),
      width = lget(sizing, "width", NA),
      height = lget(sizing, "height", NA),
      units = lget(sizing, "units", "in")
    )
  } else {
    plot
  }
}

#' Return a list mapping values to colors in `palette`
#'
#' @param palette The name of a paletteer_d palette
map_colors_d <- function(values, palette) {
  cmap <- paletteer_d(palette)
  lengths <- map_dbl(list(values, cmap), length)
  if (lengths[1] > lengths[2]) {
    msg <- glue(
      "The number of values exceeds the number of colors in {palette}!"
    )
    msg <- glue("{msg}\n N values: {lengths[1]}, N colors: {lengths[2]} ")
    stop(msg)
  }
  p <- cmap[1:lengths[1] + 1] |> as.list()
  setNames(p, values)
}


pca_dgelist <- function(
  dgelist,
  plot_aes = list(shape = "group", color = "type"),
  log = TRUE,
  to_cpm = TRUE
) {
  if (class(dgelist) == "DGEList" && to_cpm) {
    counts <- edgeR::cpm(dgelist$counts, log = log) # Convert counts into counts per million
    # When dgelist is normalized with `normLibSizes`, then cpm automatically generates
    # normalized counts from the computed normalization factors
    # This is roughly equivalent to the following:
    # 1. Get effective library size
    ## effective_lib_size <- dgelist$samples$lib.size * dgelist$samples$norm.factors
    # 2. Scale the counts by the lib size
    ## counts <- t(t(counts) / effective_lib_size)
  } else {
    counts <- dgelist$counts
  }
  pca <- prcomp(t(counts))
  summary <- summary(pca)
  var_pca1 <- summary$importance[, 1]["Proportion of Variance"]
  var_pca2 <- summary$importance[, 2]["Proportion of Variance"]
  pcs <- as_tibble(pca$x) |> bind_cols(dgelist$samples)
  plot_aes <- rlang::syms(plot_aes)
  aes_fn <- do.call(\(...) aes(x = PC1, y = PC2, ...), plot_aes)
  plot <- ggplot(pcs, aes_fn) +
    geom_point() +
    xlab(glue("PC1 ({var_pca1})")) +
    ylab(glue("PC2 ({var_pca2})"))
  plot
}


## * Gene-level variant plot

get_range_list <- function(gr, offset = TRUE) {
  if (offset) {
    gr <- sort(gr)
    o <- -start(gr[1]) + 1
    gr <- shift(gr, o)
  } else {
    o <- NULL
  }
  ranges <- lapply(seq_along(gr), \(i) as.integer(c(start(gr[i]), end(gr[i]))))
  return(list(ranges = ranges, offset = o))
}

#' Add an exon track with variant data to GenomeViz object `gv`
#'
#' @param exons GRanges object with exon ranges
#' @param variants GRanges object  containing variants (must have mcols containing
#' the variant information)
#' @param gene GRanges object with 1 range containing the gene
add_var_exon_track <- function(
  gv,
  gene,
  exons,
  variants,
  cmap,
  track_name = NULL
) {
  ranges <- get_range_list(exons)
  if (!is.null(track_name)) {
    name <- track_name
  } else {
    name <- names(gene)
  }
  track <- gv$add_feature_track(name, end(gene) - start(gene))
  rl <- get_range_list(exons, TRUE)

  strand <- ifelse(all(strand(exons) == "-"), -1, 1)
  track$add_exon_feature(rl$ranges, strand = strand, plotstyle = "box")

  for (i in seq_along(variants)) {
    v <- variants[i]
    loc <- as.integer(start(v) + rl$offset)
    var_names <- str_split_1(v$existing_variation, ";") # Prefer cosmic
    only_cosmic <- purrr::keep(var_names, \(x) str_detect(x, "COSV"))
    if (length(only_cosmic) > 0) {
      vn <- first(only_cosmic)
    } else {
      vn <- first(var_names)
    }
    color <- tryCatch(
      expr = {
        cmap[[first(v$consequence)]]
      },
      error = \(cnd) {
        print("WARNING, failed to get color for v")
        print(cmap)
        print(v)
        NULL
      }
    )
    if (is.null(color)) {
      next
    }
    label <- glue("{v$ref}>{v$alt} {vn}")
    track$add_text(
      x = loc,
      text = label,
      fontname = "monospace",
      backgroundcolor = color,
      show_line = TRUE,
      line_kws = list(color = "black", linewidth = "1")
    )
  }
}

#' Produce a plot comparing the sample
#'
#' @param ensdb Ensembldb object
#'
plot_sample_variants <- function(
  ensdb,
  vep_files,
  gene_name,
  outfile,
  canonical_tx = NULL,
  sample_names = NULL,
  allowed_consequences = c(
    "missense_variant",
    "frameshift_variant",
    "downstream_gene_variant",
    "upstream_gene_variant",
    "stop_gained",
    "splice_region_variant",
    "inframe_deletion",
    "splice_donor_5th_base_variant"
  ),
  palette = ""
) {
  library(ensembldb)
  pyplt <- new.env()
  reticulate::source_python(paste0(PY_SRC, "/", "plotting.py"), envir = pyplt)
  pg <- reticulate::import("pygenomeviz")

  transcripts <- exonsBy(ensdb, "tx", filter = GeneNameFilter(gene_name))
  gene <- genes(db, filter = ~ symbol == gene_name)
  if (length(gene) == 0) {
    stop(glue("Gene {gene_name} not found in ensdb"))
  }
  if (
    is.null(canonical_tx) || intersect(names(transcripts), canonical_tx) <= 0
  ) {
    chosen_tx <- transcripts[[1]]
  } else {
    chosen_tx <- transcripts[names(transcripts) %in% canonical_tx][[1]]
  }
  gv <- pg$GenomeViz()
  if (is.null(sample_names)) {
    sample_names <- basename_no_ext(vep_files)
  }

  vep_tbs <- lapply(vep_files, read_tsv)

  all_consequences <- bind_rows(vep_tbs)$Consequence |>
    unique() |>
    flatten_by(collapse = FALSE, unique = TRUE) |>
    unlist() |>
    keep(\(x) x %in% allowed_consequences)

  cmap <- map_colors_d(all_consequences, palette)

  seen_consequence <- c()

  plot_helper <- function(file, name) {
    gr <- into_granges(file, allowed_consequences = allowed_consequences)
    cur_vars <- gr[replace_na(gr$symbol == gene_name, FALSE)]
    seen_consequence <<- c(
      seen_consequence,
      unique(unlist(mcols(cur_vars)$consequence))
    )
    add_var_exon_track(
      gv,
      gene,
      chosen_tx,
      cur_vars,
      track_name = name,
      cmap = cmap
    )
  }

  for (i in seq_along(vep_files)) {
    plot_helper(vep_tbs[[i]], sample_names[i])
  }
  seen_consequence <- unique(seen_consequence)
  cmap <- cmap[names(cmap) %in% seen_consequence]

  pyplt$savefig(gv, outfile, custom_legend = cmap)
}

## * Graph visualizations

prepare_tflink <- function(gene_reference, mitab_file) {
  id_map <- gene_reference
  tflink <- simplify_tflink_mitab(mitab_file)
  ncbi2hgnc <- setNames(id_map$hgnc_symbol, id_map$entrezgene_id)
  tflink$regulator <- ncbi2hgnc[tflink$regulator]
  tflink$target <- ncbi2hgnc[tflink$target]
  tflink <- filter(tflink, !(is.na(target) | is.na(regulator)))
  as_tbl_graph(tflink, directed = TRUE) |>
    activate(edges) |>
    mutate(
      ends = "last",
      elabel = replace_values(
        regulator_type,
        NA ~ "",
        "activator" ~ "activates",
        "repressor" ~ "represses",
        "activator|repressor" ~ "activates|represses"
      )
    )
}

prepare_reactome_pwy <- function(pathway_file, relation_file) {
  pwys <- read_tsv(
    pathway_file,
    col_names = c("id", "name", "species")
  )
  relations <- read_tsv(relation_file, col_names = c("from", "to"))
  as_tbl_graph(relations) |>
    activate(nodes) |>
    rename(id = "name") |>
    inner_join(pwys, by = join_by(id)) |>
    filter(species == "Homo sapiens")
}

prepare_reactome_fi <- function(file) {
  read_tsv(file) |>
    filter(!str_detect(Annotation, "expression") & Annotation != "predicted") |>
    as_tbl_graph() |>
    activate(edges) |>
    mutate(
      ends = replace_values(
        Direction,
        c("<-", "|-", "<-|") ~ "first",
        c("->", "-|", "|->") ~ "last",
        c("<->", "|-|") ~ "both",
        "-" ~ "last" # This is wrong, but we can remove the arrow by setting angle to 0
      )
    ) |>
    rename(elabel = "Annotation")
}

plot_de_graph <- function(G, palette_c = "grDevices::Cyan-Magenta") {
  ends <- E(G)$ends
  activate(G, nodes) |>
    mutate(lfc = replace_values(lfc, !is_de ~ NA)) |>
    ggraph("kk") +
    geom_edge_link(
      aes(label = elabel),
      angle_calc = "along",
      label_dodge = unit(3, "mm"),
      arrow = grid::arrow(
        ends = ends,
        type = "closed",
        angle = ifelse(is.na(ends), 0, 25),
        length = unit(0.1, "inches")
      ),
      end_cap = circle(1, "cm"),
      start_cap = circle(1, "cm")
    ) +
    geom_node_label(
      aes(label = name, fill = lfc, fontface = ifelse(is_de, "bold", "plain")),
    ) +
    theme(panel.background = element_rect(fill = "white")) +
    scale_fill_paletteer_c(
      palette_c,
      na.value = "white"
    ) +
    guides(
      fill = guide_legend("LFC"),
      color = guide_legend("LFC")
    )
}

plot_enriched_graph <- function(
  G,
  palette_c = "ggthemes::Classic Red-Green Light",
  palette_d = "LaCroixColoR::CeriseLimon",
  layout = "kk"
) {
  expect_col_exists(as_tibble(G, active = "nodes"), c(stat, name, enriched))
  ggraph(G, layout) +
    geom_edge_diagonal() +
    geom_node_point(aes(color = enriched), size = 2) +
    geom_node_label(aes(label = name, fill = stat), repel = TRUE) +
    theme(panel.background = element_rect(fill = "white")) +
    scale_color_paletteer_d(palette_d) +
    scale_fill_paletteer_c(palette_c, na.value = "white") +
    guides(
      fill = guide_legend("Enrichment test statistic"),
      color = "none"
    )
}

ggsave_graph_dynamic <- function(G, plot, filename) {
  assert(
    check_class(G, "igraph"),
    check_class(plot, "ggplot"),
    check_character(filename),
    combine = "and"
  )
  if (length(G) >= 130) {
    dim <- 30
  } else if (length(G) >= 100) {
    dim <- 25
  } else if (length(G) >= 30) {
    dim <- 12
  } else if (length(G) >= 10) {
    dim <- 10
  } else if (length(G) >= 5) {
    dim <- 8
  } else {
    dim <- 5
  }
  ggsave(
    filename = filename,
    plot = plot,
    height = dim,
    width = dim,
    units = "in"
  )
}

combine_sample_label_plots <- function(
  label_tbs,
  palettes = list()
) {
  assert_list(label_tbs, names = "named")
  if (length(palettes) > 1) {
    assert_list(palettes, names = "named")
    assert_names(names(palettes), subset.of = names(label_tbs))
  }

  plot_one <- function(tb, title, lpos, last = FALSE, palette = NULL) {
    assert_data_frame(tb)
    plot <- ggplot(tb, aes(x = sample, y = 1, fill = label)) +
      TILE_CALL +
      theme_minimal() +
      theme(
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        axis.text.y = element_blank()
      )
    if (!is.null(palette)) {
      plot <- plot + scale_fill_paletteer_d(palette)
    }
    if (!last) {
      plot <- plot +
        theme(axis.title.x = element_blank(), axis.text.x = element_blank())
    }
    plot + guides(fill = guide_legend(title = title, position = lpos))
  }

  plots <- lapply(
    seq_along(label_tbs),
    \(i) {
      name <- names(label_tbs)[i]
      position <- ifelse(i %% 2 == 0, "left", "right")
      plot_one(
        label_tbs[[name]],
        name,
        palettes[[name]],
        lpos = position,
        last = i == length(label_tbs)
      ) +
        theme(plot.margin = margin(0, 0, 0, 0, "cm"))
    }
  )
  wrap_plots(plots, ncol = 1)
}
