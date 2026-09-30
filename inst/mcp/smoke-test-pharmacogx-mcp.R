#!/usr/bin/env Rscript

required_packages <- c("jsonlite", "processx")
missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]
if (length(missing_packages) > 0) {
  stop(
    "MCP transport smoke testing requires: ",
    paste(missing_packages, collapse = ", ")
  )
}

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_args, value = TRUE)
script_dir <- dirname(normalizePath(
  sub("^--file=", "", file_arg[[1]]),
  mustWork = TRUE
))
runner <- file.path(script_dir, "run-pharmacogx-mcp.R")

process <- processx::process$new(
  Sys.which("Rscript"),
  runner,
  env = c(
    PHARMACOGX_MCP_USE_SOURCE_PACKAGE = "true",
    PHARMACOGX_MCP_SINGLE_CONTENT = "true",
    PHARMACOGX_MCP_SESSION_TOOLS = "false"
  ),
  stdin = "|",
  stdout = "|",
  stderr = "|",
  cleanup = TRUE
)
on.exit(process$kill(), add = TRUE)

send_message <- function(message) {
  process$write_input(paste0(
    jsonlite::toJSON(message, auto_unbox = TRUE, null = "null"),
    "\n"
  ))
}

read_response <- function(id, timeout_ms = 30000) {
  started <- Sys.time()
  repeat {
    process$poll_io(250)
    lines <- process$read_output_lines()
    for (line in lines) {
      parsed <- tryCatch(
        jsonlite::fromJSON(line, simplifyVector = FALSE),
        error = function(e) NULL
      )
      if (!is.null(parsed) && identical(parsed$id, id)) {
        return(parsed)
      }
    }
    if (!process$is_alive()) {
      stop(
        "MCP server exited before response ",
        id,
        ": ",
        paste(process$read_error_lines(), collapse = "\n")
      )
    }
    elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs")) * 1000
    if (elapsed > timeout_ms) {
      stop("Timed out waiting for MCP response ", id, ".")
    }
  }
}

send_message(list(
  jsonrpc = "2.0",
  id = 1L,
  method = "initialize",
  params = list(
    protocolVersion = "2025-11-25",
    capabilities = list(),
    clientInfo = list(name = "pharmacogx-smoke", version = "1.0")
  )
))
initialize <- read_response(1L)
if (!identical(initialize$result$protocolVersion, "2025-11-25")) {
  stop("Unexpected MCP protocol version.")
}

send_message(list(
  jsonrpc = "2.0",
  method = "notifications/initialized",
  params = list()
))
send_message(list(
  jsonrpc = "2.0",
  id = 2L,
  method = "tools/list",
  params = list()
))
tools_response <- read_response(2L)
tool_names <- vapply(
  tools_response$result$tools,
  function(tool) tool$name,
  character(1)
)
if (length(tool_names) != 33L) {
  stop("Expected 33 PharmacoGx tools, found ", length(tool_names), ".")
}
if (any(c("list_r_sessions", "select_r_session") %in% tool_names)) {
  stop("General R session tools should be disabled by default.")
}

send_message(list(
  jsonrpc = "2.0",
  id = 3L,
  method = "tools/call",
  params = list(name = "pgx_session_info", arguments = list())
))
session_response <- read_response(3L)
if (length(session_response$result$content) != 1L) {
  stop("Expected exactly one MCP content block.")
}
if (!is.null(session_response$result$structuredContent)) {
  stop("Single-content mode should not return a second structured payload.")
}
session <- jsonlite::fromJSON(session_response$result$content[[1]]$text)
if (!identical(session$tool_count, 33L)) {
  stop("Session health response reported the wrong tool count.")
}

message(
  "PharmacoGx MCP transport smoke test passed: ",
  length(tool_names),
  " tools, one JSON payload, source=",
  session$pharmacogx_package_source
)
