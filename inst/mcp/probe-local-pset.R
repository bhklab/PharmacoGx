#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop(
    "Usage: probe-local-pset.R <tool-file> <package-source> ",
    "<pset-path> <result-path>"
  )
}

tool_file <- normalizePath(args[[1]], mustWork = TRUE)
package_source <- args[[2]]
pset_path <- normalizePath(args[[3]], mustWork = TRUE)
result_path <- args[[4]]

result <- tryCatch(
  {
    use_source <- !identical(package_source, "installed") &&
      file.exists(file.path(package_source, "DESCRIPTION"))
    if (use_source) {
      if (!requireNamespace("pkgload", quietly = TRUE)) {
        stop("Isolated source-checkout probing requires the pkgload package.")
      }
      suppressMessages(suppressWarnings(pkgload::load_all(
        package_source,
        attach = FALSE,
        helpers = FALSE,
        quiet = TRUE,
        warn_conflicts = FALSE,
        debug = FALSE
      )))
      Sys.setenv(PHARMACOGX_MCP_PACKAGE_SOURCE = package_source)
    } else {
      if (!requireNamespace("PharmacoGx", quietly = TRUE)) {
        stop("The PharmacoGx package is unavailable in the isolated worker.")
      }
      Sys.setenv(PHARMACOGX_MCP_PACKAGE_SOURCE = "installed")
    }

    Sys.setenv(
      PHARMACOGX_MCP_IN_PROBE_WORKER = "true",
      PHARMACOGX_MCP_TOOL_FILE = tool_file,
      PHARMACOGX_MCP_DIR = dirname(tool_file)
    )
    tool_environment <- new.env(parent = globalenv())
    source(tool_file, local = tool_environment)
    snapshot <- tool_environment$pgx_build_dataset_capability_snapshot(
      pset_path
    )
    list(ok = TRUE, snapshot = snapshot)
  },
  error = function(e) {
    list(ok = FALSE, message = conditionMessage(e))
  }
)

saveRDS(result, result_path, compress = "gzip")
if (!isTRUE(result$ok)) {
  message(result$message)
  quit(status = 1L)
}
