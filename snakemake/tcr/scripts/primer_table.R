library(gt)
library(glue)
library(tidyverse)

make_chain_tb <- function(lst, chain) {
  cur <- lst[[chain]]
  table <- lapply(c("pair 1", "pair 2"), \(name) {
    pair_info <- cur[[name]]
    tb <- lapply(
      c("Forward", "Reverse"),
      \(x) {
        tibble(
          direction = x,
          gc = pair_info[[x]]$gc_content,
          length = pair_info[[x]]$length,
          seq = pair_info[[x]]$seq,
          tm = pair_info[[x]]$tm,
          pair = name
        )
      }
    ) |>
      bind_rows()
    tb |>
      add_row(
        seq = pair_info$figure,
        .after = 2,
        pair = name
      )
  }) |>
    bind_rows()
}


tab_formatting <- function(tab) {
  tab |>
    fmt_passthrough("seq") |>
    fmt_number(columns = c("gc", "tm")) |>
    cols_hide("chain") |>
    sub_missing(missing_text = "") |>
    cols_align("seq", align = "center") |>
    tab_style(
      style = cell_text(font = "Monospace"),
      locations = cells_body(columns = "seq")
    ) |>
    cols_label_with(fn = \(x) str_to_title(x)) |>
    cols_label(seq = "Sequence", gc = "GC content", tm = "TM (celsius)") |>
    tab_style(
      style = cell_text(
        font = "Monospace",
        whitespace = "break-spaces",
        align = "left"
      ),
      locations = cells_body(columns = "seq", rows = is.na(length))
    )
}

make_primers_table <- function(chain_tb, chain_data) {
  subtitle_tmp <- c()
  subtitle <- map_chr(c("v", "d", "j"), \(gene) {
    val <- chain_data[[glue("{gene}_call")]]
    if (!is.null(val)) {
      glue("<b>{str_to_upper(gene)}</b> call: {val}")
    } else {
      ""
    }
  }) |>
    discard(\(x) x == "") |>
    paste0(collapse = ", ")
  gt(chain_tb, rowname_col = "direction", groupname_col = "pair") |>
    tab_formatting() |>
    tab_header(title = "", subtitle = html(subtitle))
}


make_sample_primer_table <- function(sample_data, sample_key = "Sample_Name") {
  combined <- lapply(
    c("VJ_1", "VDJ_1"),
    \(c) make_chain_tb(sample_data, c) |> mutate(chain = c)
  ) |>
    bind_rows() |>
    mutate(pair = str_to_title(pair))
  sname <- sample_data[[sample_key]]
  title <- glue("{sname} fusion primers")
  subtitle <- glue(
    "Clone ID: {sample_data$clone_id}, Chain pairing: {sample_data$chain_pairing}"
  )
  gt(
    combined,
    rowname_col = c("pair", "direction"),
    groupname_col = "chain",
    process_md = TRUE
  ) |>
    tab_formatting() |>
    opt_stylize(style = 6, add_row_striping = FALSE, color = "green") |>
    text_transform(
      \(clist) {
        chain <- head(unlist(clist), n = 1)
        chain_data <- sample_data[[chain]]
        gene_info <- map_chr(c("v", "d", "j"), \(gene) {
          val <- chain_data[[glue("{gene}_call")]]
          if (!is.null(val)) {
            glue("<b>{str_to_upper(gene)}</b> call: {val}")
          } else {
            ""
          }
        }) |>
          discard(\(x) x == "") |>
          paste0(collapse = ", ")
        html(glue("{chain}: {gene_info}"))
      },
      locations = cells_row_groups()
    ) |>
    text_transform(
      locations = cells_row_groups(),
      fn = \(x) lapply(x, html)
    ) |>
    tab_header(title = title, subtitle = subtitle)
}

write_sample_primers_fasta <- function(sample_data, outfile) {
  sname <- sample_data$Sample_Name
  clone_id <- sample_data$clone_id
  entries <- lapply(c("VJ_1", "VDJ_1"), \(chain) {
    chain_data <- sample_data[[chain]]
    lapply(c("pair 1", "pair 2"), \(pair) {
      pname <- str_replace(pair, "pair ", "p")
      directions <- c("Forward", "Reverse")
      seqs <- map_chr(directions, \(d) chain_data[[pair]][[d]]$seq)
      glue(
        ">{sname}_cid{clone_id}_{chain}_{pname}_{str_to_lower(directions)}\n{seqs}"
      )
    })
  }) |>
    unlist(use.names = FALSE)
  write_lines(entries, outfile)
}


write_primer_tables <- function() {
  data <- jsonlite::read_json(snakemake@input[[1]])
  outdir <- snakemake@params[["outdir"]]
  dir.create(outdir)
  for (d in data) {
    outfile <- glue("{outdir}/cid{d$clone_id}_{d$Sample_Name}.html")
    fasta_out <- glue("{outdir}/cid{d$clone_id}_{d$Sample_Name}.fasta")
    table <- make_sample_primer_table(d)
    write_sample_primers_fasta(d, fasta_out)
    as_raw_html(table, inline_css = TRUE) |> htmltools::save_html(outfile)
  }
}

if (exists("snakemake")) {
  globalenv()[[snakemake@rule]]()
}
