#!/usr/bin/env Rscript

if (!requireNamespace("yaml", quietly = TRUE)) {
  stop("Case validation requires the yaml package.")
}
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("Case validation requires the jsonlite package.")
}

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_args, value = TRUE)
script_dir <- if (length(file_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE))
} else {
  getwd()
}
default_manifest <- file.path(dirname(script_dir), "development-cases.yaml")
args <- commandArgs(trailingOnly = TRUE)
manifest_path <- if (length(args) > 0) args[[1]] else default_manifest
manifest_path <- normalizePath(manifest_path, mustWork = TRUE)
suite <- yaml::read_yaml(manifest_path)

errors <- character()
required_suite <- c("schema_version", "suite_id", "case_set", "cases")
missing_suite <- setdiff(required_suite, names(suite))
if (length(missing_suite) > 0) {
  errors <- c(
    errors,
    paste("Missing suite fields:", paste(missing_suite, collapse = ", "))
  )
}
if (!identical(suite$schema_version, 1L)) {
  errors <- c(errors, "schema_version must be 1")
}
if (!suite$case_set %in% c("development", "held_out")) {
  errors <- c(errors, "case_set must be development or held_out")
}
if (!is.list(suite$cases) || length(suite$cases) == 0) {
  errors <- c(errors, "cases must be a non-empty list")
}

required_case <- c(
  "case_id",
  "task_class",
  "prompt",
  "data_scope",
  "reference",
  "scoring",
  "limits"
)
allowed_statuses <- c("draft", "needs_independent_review", "verified")
case_ids <- character()
for (i in seq_along(suite$cases)) {
  case <- suite$cases[[i]]
  label <- paste0("case[", i, "]")
  missing_case <- setdiff(required_case, names(case))
  if (length(missing_case) > 0) {
    errors <- c(
      errors,
      paste(label, "is missing", paste(missing_case, collapse = ", "))
    )
    next
  }
  case_ids <- c(case_ids, case$case_id)
  if (!grepl("^[a-z0-9-]+$", case$case_id)) {
    errors <- c(errors, paste(label, "has an invalid case_id"))
  }
  if (!nzchar(trimws(case$prompt))) {
    errors <- c(errors, paste(label, "has an empty prompt"))
  }
  if (!case$reference$status %in% allowed_statuses) {
    errors <- c(errors, paste(label, "has an invalid reference status"))
  }
  if (length(case$reference$hard_failures) == 0) {
    errors <- c(errors, paste(label, "must define at least one hard failure"))
  }
  if (length(case$scoring$critical_dimensions) == 0) {
    errors <- c(errors, paste(label, "must define critical dimensions"))
  }
  if (case$limits$max_tool_calls < 1 || case$limits$max_runtime_seconds < 1) {
    errors <- c(errors, paste(label, "has invalid execution limits"))
  }
}
duplicates <- unique(case_ids[duplicated(case_ids)])
if (length(duplicates) > 0) {
  errors <- c(
    errors,
    paste("Duplicate case IDs:", paste(duplicates, collapse = ", "))
  )
}

summary <- list(
  valid = length(errors) == 0,
  manifest = manifest_path,
  suite_id = suite$suite_id,
  case_set = suite$case_set,
  case_count = length(suite$cases),
  verified_references = sum(vapply(
    suite$cases,
    function(case) identical(case$reference$status, "verified"),
    logical(1)
  )),
  errors = errors
)
cat(jsonlite::toJSON(summary, auto_unbox = TRUE, pretty = TRUE), "\n")
if (length(errors) > 0) {
  quit(status = 1L)
}
