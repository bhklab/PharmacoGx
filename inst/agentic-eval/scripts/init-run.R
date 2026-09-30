#!/usr/bin/env Rscript

if (!requireNamespace("yaml", quietly = TRUE)) {
  stop("Run initialization requires the yaml package.")
}
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("Run initialization requires the jsonlite package.")
}

parse_args <- function(args) {
  result <- list()
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--") || i == length(args)) {
      stop("Arguments must be provided as --name value pairs.")
    }
    result[[sub("^--", "", key)]] <- args[[i + 1L]]
    i <- i + 2L
  }
  result
}

split_values <- function(value) {
  if (is.null(value) || !nzchar(value)) {
    return(character())
  }
  values <- trimws(strsplit(value, ",", fixed = TRUE)[[1]])
  values[nzchar(values)]
}

options <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c("host", "condition", "model", "plugin-commit", "output")
missing <- required[
  !vapply(
    required,
    function(name) {
      !is.null(options[[name]]) && nzchar(options[[name]])
    },
    logical(1)
  )
]
if (length(missing) > 0) {
  stop("Missing required arguments: ", paste(missing, collapse = ", "))
}
allowed_conditions <- c("baseline", "skill", "mcp", "plugin")
if (!options$condition %in% allowed_conditions) {
  stop(
    "--condition must be one of: ",
    paste(allowed_conditions, collapse = ", "),
    "."
  )
}

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_args, value = TRUE)
script_dir <- if (length(file_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE))
} else {
  getwd()
}
default_manifest <- file.path(dirname(script_dir), "development-cases.yaml")
default_skill <- file.path(
  dirname(dirname(script_dir)),
  "skills",
  "pharmacogx-agentic-analysis",
  "SKILL.md"
)
manifest_path <- if (is.null(options$manifest)) {
  default_manifest
} else {
  options$manifest
}
manifest_path <- normalizePath(manifest_path, mustWork = TRUE)
suite <- yaml::read_yaml(manifest_path)
uses_skill <- options$condition %in% c("skill", "plugin")
skill_path <- if (!uses_skill) {
  NULL
} else if (is.null(options[["skill-path"]])) {
  default_skill
} else {
  options[["skill-path"]]
}
skill <- if (is.null(skill_path)) {
  NULL
} else {
  normalized_skill_path <- normalizePath(skill_path, mustWork = TRUE)
  list(
    name = "pharmacogx-agentic-analysis",
    path = normalized_skill_path,
    md5 = unname(tools::md5sum(normalized_skill_path))
  )
}
resource_paths <- split_values(options[["resource-paths"]])
attached_resources <- lapply(resource_paths, function(path) {
  normalized_path <- normalizePath(path, mustWork = TRUE)
  list(
    name = basename(normalized_path),
    path = normalized_path,
    md5 = unname(tools::md5sum(normalized_path))
  )
})
output_path <- normalizePath(options$output, mustWork = FALSE)
if (file.exists(output_path) && !identical(options$overwrite, "true")) {
  stop("Output already exists. Pass --overwrite true to replace it.")
}

created_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
execution_routes <- c(
  baseline = "host_native",
  skill = "skill_direct",
  mcp = "mcp",
  plugin = "skill_mcp"
)
run_id <- paste(
  suite$suite_id,
  options$host,
  options$condition,
  format(Sys.time(), "%Y%m%dT%H%M%S%z"),
  sep = "-"
)
run <- list(
  schema_version = 1L,
  run_id = run_id,
  created_at = created_at,
  suite = list(
    suite_id = suite$suite_id,
    case_set = suite$case_set,
    manifest = manifest_path,
    manifest_md5 = unname(tools::md5sum(manifest_path))
  ),
  environment = list(
    host = options$host,
    condition = options$condition,
    model = options$model,
    model_version = options[["model-version"]],
    plugin_commit = options[["plugin-commit"]],
    execution_route = unname(execution_routes[[options$condition]]),
    mcp_enabled = options$condition %in% c("mcp", "plugin"),
    skill = skill,
    attached_resources = attached_resources,
    temperature = if (is.null(options$temperature)) {
      NULL
    } else {
      as.numeric(options$temperature)
    },
    pset_versions = structure(list(), names = character())
  ),
  cases = lapply(suite$cases, function(case) {
    list(
      case_id = case$case_id,
      prompt = case$prompt,
      status = "pending",
      raw_response = NULL,
      tool_trace = list(),
      artifacts = list(),
      timing = list(started_at = NULL, completed_at = NULL, seconds = NULL),
      adjudication = list(
        reviewer = NULL,
        dimension_scores = list(),
        automatic_failures = list(),
        passed = NULL,
        notes = NULL
      )
    )
  })
)

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
jsonlite::write_json(
  run,
  output_path,
  auto_unbox = TRUE,
  pretty = TRUE,
  null = "null"
)
cat(output_path, "\n")
