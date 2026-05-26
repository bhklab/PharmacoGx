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

tool_file <- system.file("mcp", "PharmacoGx-tools.R", package = "PharmacoGx")
if (!nzchar(tool_file)) {
  tool_file <- file.path(script_dir, "PharmacoGx-tools.R")
}
if (!file.exists(tool_file)) {
  stop("Could not find PharmacoGx-tools.R next to the MCP runner script.")
}

transport <- Sys.getenv("PHARMACOGX_MCP_TRANSPORT", "stdio")
port <- as.integer(Sys.getenv("MCPTOOLS_PORT", "8080"))

mcptools::mcp_server(
  tools = tool_file,
  type = transport,
  port = port
)
