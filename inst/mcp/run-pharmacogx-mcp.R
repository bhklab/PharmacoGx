#!/usr/bin/env Rscript

if (!requireNamespace("mcptools", quietly = TRUE)) {
  stop(
    "The PharmacoGx MCP demo requires mcptools. Install it with ",
    "pak::pkg_install('mcptools')."
  )
}

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_args, value = TRUE)
script_dir <- if (length(file_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE))
} else {
  getwd()
}

source_tool_file <- file.path(script_dir, "PharmacoGx-tools.R")
installed_tool_file <- system.file(
  "mcp",
  "PharmacoGx-tools.R",
  package = "PharmacoGx"
)
tool_file <- if (file.exists(source_tool_file)) {
  source_tool_file
} else {
  installed_tool_file
}
if (!file.exists(tool_file)) {
  stop("Could not find PharmacoGx-tools.R next to the MCP runner script.")
}

package_root <- normalizePath(
  file.path(script_dir, "..", ".."),
  mustWork = FALSE
)
use_source_package <- tolower(Sys.getenv(
  "PHARMACOGX_MCP_USE_SOURCE_PACKAGE",
  "auto"
))
source_package_available <- file.exists(file.path(package_root, "DESCRIPTION"))
should_load_source <- identical(use_source_package, "true") ||
  (identical(use_source_package, "auto") && source_package_available)
if (should_load_source) {
  if (!source_package_available) {
    stop(
      "PHARMACOGX_MCP_USE_SOURCE_PACKAGE=true, but no DESCRIPTION was found at ",
      package_root
    )
  }
  if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop(
      "Loading PharmacoGx from this source checkout requires pkgload. Install ",
      "it with pak::pkg_install('pkgload')."
    )
  }
  suppressMessages(suppressWarnings(pkgload::load_all(
    package_root,
    attach = FALSE,
    helpers = FALSE,
    quiet = TRUE,
    warn_conflicts = FALSE,
    debug = FALSE
  )))
  Sys.setenv(PHARMACOGX_MCP_PACKAGE_SOURCE = package_root)
} else {
  Sys.setenv(PHARMACOGX_MCP_PACKAGE_SOURCE = "installed")
}

Sys.setenv(PHARMACOGX_MCP_DIR = dirname(tool_file))
Sys.setenv(
  PHARMACOGX_MCP_TOOL_FILE = normalizePath(tool_file, mustWork = TRUE),
  PHARMACOGX_MCP_SERVER_STARTED_AT = format(
    Sys.time(),
    "%Y-%m-%dT%H:%M:%S%z"
  )
)

log_file <- Sys.getenv("PHARMACOGX_MCP_LOG_FILE", "")
if (!nzchar(log_file)) {
  log_file <- file.path(
    tools::R_user_dir("PharmacoGx", which = "cache"),
    "mcp.log"
  )
  Sys.setenv(PHARMACOGX_MCP_LOG_FILE = log_file)
}
dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)
write(
  paste(
    format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    "server_start",
    paste0("pid=", Sys.getpid()),
    paste0("transport=", Sys.getenv("PHARMACOGX_MCP_TRANSPORT", "stdio")),
    paste0("package_source=", Sys.getenv("PHARMACOGX_MCP_PACKAGE_SOURCE")),
    paste0("tool_file=", normalizePath(tool_file, mustWork = TRUE)),
    sep = "\t"
  ),
  file = log_file,
  append = TRUE
)

transport <- Sys.getenv("PHARMACOGX_MCP_TRANSPORT", "stdio")
port <- as.integer(Sys.getenv("MCPTOOLS_PORT", "8080"))

mcptools::mcp_server(
  tools = tool_file,
  type = transport,
  port = port
)
