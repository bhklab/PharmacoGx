#!/usr/bin/env Rscript

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("Snapshot export requires the jsonlite package.")
}
if (!requireNamespace("ellmer", quietly = TRUE)) {
  stop("Snapshot export requires the ellmer package.")
}

parse_args <- function(args) {
  parsed <- list()
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--") || i == length(args)) {
      stop("Arguments must be provided as --name value pairs.")
    }
    parsed[[sub("^--", "", key)]] <- args[[i + 1L]]
    i <- i + 2L
  }
  parsed
}

split_values <- function(value) {
  if (is.null(value) || !nzchar(value)) {
    return(NULL)
  }
  values <- trimws(strsplit(value, ",", fixed = TRUE)[[1]])
  values[nzchar(values)]
}

parse_flag <- function(value, default = FALSE) {
  if (is.null(value)) {
    return(default)
  }
  tolower(value) %in% c("1", "true", "yes", "on")
}

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_args, value = TRUE)
script_dir <- if (length(file_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE))
} else {
  getwd()
}
default_package_source <- dirname(dirname(dirname(script_dir)))

cli <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c("datasets", "output")
missing <- required[
  !vapply(
    required,
    function(name) !is.null(cli[[name]]) && nzchar(cli[[name]]),
    logical(1)
  )
]
if (length(missing) > 0) {
  stop("Missing required arguments: ", paste(missing, collapse = ", "))
}

package_source <- normalizePath(
  if (is.null(cli[["package-source"]])) {
    default_package_source
  } else {
    cli[["package-source"]]
  },
  mustWork = TRUE
)
tool_file <- file.path(package_source, "inst", "mcp", "PharmacoGx-tools.R")
tool_file <- normalizePath(tool_file, mustWork = TRUE)
output_path <- normalizePath(cli$output, mustWork = FALSE)
if (file.exists(output_path) && !parse_flag(cli$overwrite)) {
  stop("Output already exists. Pass --overwrite true to replace it.")
}

if (!is.null(cli[["local-pset-dirs"]])) {
  Sys.setenv(PHARMACOGX_LOCAL_PSET_DIRS = cli[["local-pset-dirs"]])
}
if (!is.null(cli[["local-pset-manifest"]])) {
  Sys.setenv(PHARMACOGX_LOCAL_PSET_MANIFEST = cli[["local-pset-manifest"]])
}
Sys.setenv(
  PHARMACOGX_MCP_PACKAGE_SOURCE = package_source,
  PHARMACOGX_MCP_TOOL_FILE = tool_file,
  PHARMACOGX_MCP_DIR = dirname(tool_file)
)

tool_environment <- new.env(parent = globalenv())
source(tool_file, local = tool_environment)

datasets <- split_values(cli$datasets)
result <- tool_environment$pgx_check_analysis_requirements(
  datasets = datasets,
  drug = cli$drug,
  sample_query = cli[["sample-query"]],
  metric = cli$metric,
  mDataType = cli[["mdata-type"]],
  features = split_values(cli$features),
  include_local = TRUE,
  include_bundled = parse_flag(cli[["include-bundled"]]),
  limit = length(datasets)
)

compute_md5 <- parse_flag(cli[["compute-md5"]])
files <- lapply(datasets, function(dataset) {
  path <- tool_environment$pgx_local_dataset_path(dataset)
  if (is.null(path)) {
    return(list(dataset = dataset, available = FALSE))
  }
  info <- file.info(path)
  list(
    dataset = dataset,
    available = TRUE,
    file = basename(path),
    bytes = unname(info$size),
    modified_at = format(info$mtime, "%Y-%m-%dT%H:%M:%S%z"),
    md5 = if (compute_md5) unname(tools::md5sum(path)) else NULL
  )
})

git_value <- function(args) {
  tryCatch(
    system2(
      "git",
      c("-C", package_source, args),
      stdout = TRUE,
      stderr = FALSE
    ),
    error = function(e) character()
  )
}
git_commit <- git_value(c("rev-parse", "HEAD"))
git_status <- git_value(c("status", "--short"))

snapshot <- list(
  schema_version = 1L,
  generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  reference_status = "needs_independent_review",
  source = list(
    git_commit = if (length(git_commit) > 0) git_commit[[1]] else NULL,
    worktree_dirty = length(git_status) > 0,
    tool_file = basename(tool_file),
    tool_file_md5 = unname(tools::md5sum(tool_file)),
    check_mode = "isolated_worker"
  ),
  files = files,
  result = result,
  interpretation_limits = c(
    "This is a dated capability snapshot, not a biomarker result.",
    "Unchecked datasets must not be reported as unavailable.",
    "Reuse requires matching PSet file identity.",
    "Capability sample counts precede complete-case filtering."
  )
)

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
jsonlite::write_json(
  snapshot,
  output_path,
  auto_unbox = TRUE,
  pretty = TRUE,
  null = "null",
  dataframe = "rows",
  digits = NA
)
cat(output_path, "\n")
