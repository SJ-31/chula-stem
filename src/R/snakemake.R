Snakemake <- setClass(
  "Snakemake",
  slots = c(
    input = "list",
    output = "list",
    params = "list",
    wildcards = "list",
    threads = "numeric",
    log = "list",
    resources = "list",
    config = "list",
    rule = "character",
    bench_iteration = "numeric",
    scriptdir = "character"
  ),
  prototype = methods::prototype(
    input = list(),
    output = list(),
    params = list(),
    wildcards = list(),
    threads = NA_integer_,
    log = list(),
    resources = list(),
    config = list(),
    rule = NA_character_,
    bench_iteration = NA_integer_,
    scriptdir = NA_character_
  )
)
