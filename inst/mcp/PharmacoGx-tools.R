# PharmacoGx MCP tool definitions.
#
# This file is sourced by mcptools::mcp_server(tools = ...). It must evaluate to
# a named list of ellmer::tool() objects.

.datatable.aware <- TRUE

pgx_supported_datasets <- c("GDSCsmall", "CCLEsmall", "CMAPsmall")
pgx_dataset_cache <- new.env(parent = emptyenv())
pgx_capability_cache <- new.env(parent = emptyenv())

pgx_env_flag <- function(name, default = FALSE) {
  value <- Sys.getenv(name, if (isTRUE(default)) "true" else "false")
  tolower(value) %in% c("1", "true", "yes", "on")
}

pgx_normalize_text <- function(x) {
  x <- as.character(x)
  normalized <- iconv(x, from = "", to = "UTF-8", sub = "")
  normalized[is.na(normalized) & !is.na(x)] <- ""
  normalized
}

pgx_casefold <- function(x) {
  tolower(pgx_normalize_text(x))
}

pgx_error_result <- function(tool_name, error) {
  message <- conditionMessage(error)
  code <- if (grepl("Unsupported dataset|dataset must be", message)) {
    "DATASET_NOT_FOUND"
  } else if (grepl("Unknown .*field|missing required column", message)) {
    "FIELD_NOT_FOUND"
  } else if (grepl("Molecular profile|requested features", message)) {
    "MOLECULAR_DATA_UNAVAILABLE"
  } else if (grepl("Sensitivity measure", message)) {
    "METRIC_UNAVAILABLE"
  } else if (grepl("assay", message, ignore.case = TRUE)) {
    "ASSAY_UNAVAILABLE"
  } else {
    "PGX_TOOL_ERROR"
  }
  list(
    ok = FALSE,
    available = FALSE,
    error = list(
      code = code,
      tool = tool_name,
      message = message,
      remediation = paste(
        "Call pgx_session_info() to verify the active server, then use the",
        "relevant discovery or requirements-preflight tool before retrying."
      )
    )
  )
}

pgx_as_mcp_result <- function(value) {
  if (!pgx_env_flag("PHARMACOGX_MCP_SINGLE_CONTENT", default = TRUE)) {
    return(value)
  }

  json <- jsonlite::toJSON(
    value,
    auto_unbox = TRUE,
    dataframe = "rows",
    matrix = "rowmajor",
    null = "null",
    na = "null",
    digits = NA,
    pretty = FALSE
  )
  ellmer::ContentToolResult(
    value = list(ellmer::ContentText(as.character(json))),
    error = NULL,
    extra = list(),
    request = NULL
  )
}

pgx_log_file <- function() {
  log_file <- Sys.getenv("PHARMACOGX_MCP_LOG_FILE", "")
  if (nzchar(log_file)) {
    return(path.expand(log_file))
  }

  file.path(tools::R_user_dir("PharmacoGx", which = "cache"), "mcp.log")
}

pgx_log_event <- function(tool, status, detail = "") {
  log_file <- pgx_log_file()
  dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)
  detail <- gsub("[\r\n\t]+", " ", as.character(detail))
  if (nchar(detail) > 12000L) {
    detail <- paste0(substr(detail, 1L, 12000L), " [truncated]")
  }
  write(
    paste(
      format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
      tool,
      status,
      paste0("pid=", Sys.getpid()),
      detail,
      sep = "\t"
    ),
    file = log_file,
    append = TRUE
  )
  invisible(NULL)
}

pgx_execute_logged <- function(fun, tool_name, call_args, call_environment) {
  argument_names <- names(call_args)
  argument_names <- argument_names[nzchar(argument_names)]
  pgx_log_event(
    tool_name,
    "start",
    paste0("arguments=", paste(argument_names, collapse = ","))
  )

  error_calls <- NULL
  tryCatch(
    {
      value <- withCallingHandlers(
        do.call(fun, call_args, envir = call_environment),
        error = function(e) {
          error_calls <<- sys.calls()
        }
      )
      pgx_log_event(tool_name, "success")
      value
    },
    error = function(e) {
      calls <- error_calls
      if (is.null(calls)) {
        calls <- sys.calls()
      }
      calls <- vapply(
        calls,
        function(x) {
          call_text <- paste(deparse(x, width.cutoff = 200L), collapse = " ")
          if (nchar(call_text) > 500L) {
            call_text <- paste0(substr(call_text, 1L, 500L), "...")
          }
          call_text
        },
        character(1)
      )
      pgx_log_event(
        tool_name,
        "error",
        paste0(
          "message=",
          conditionMessage(e),
          " | call=",
          paste(deparse(conditionCall(e)), collapse = " "),
          " | stack=",
          paste(calls, collapse = " > ")
        )
      )
      pgx_error_result(tool_name, e)
    }
  )
}

pgx_logged_tool_function <- function(fun, tool_name) {
  wrapper <- function() NULL
  formals(wrapper) <- formals(fun)
  environment(wrapper) <- list2env(
    list(.pgx_fun = fun, .pgx_tool_name = tool_name),
    parent = environment()
  )
  body(wrapper) <- quote({
    .pgx_call_args <- as.list(match.call(expand.dots = FALSE))[-1]
    pgx_as_mcp_result(pgx_execute_logged(
      .pgx_fun,
      .pgx_tool_name,
      .pgx_call_args,
      environment()
    ))
  })
  wrapper
}

pgx_tool <- function(fun, name, description, arguments = list()) {
  writes_files <- grepl("^pgx_(download|plot)", name)
  open_world <- grepl("^pgx_(list_available_psets|download_pset)$", name)
  ellmer::tool(
    fun = pgx_logged_tool_function(fun, name),
    name = name,
    description = description,
    arguments = arguments,
    annotations = ellmer::tool_annotations(
      title = gsub("_", " ", sub("^pgx_", "PharmacoGx: ", name)),
      read_only_hint = !writes_files,
      open_world_hint = open_world,
      idempotent_hint = !identical(name, "pgx_download_pset"),
      destructive_hint = FALSE
    )
  )
}

pgx_default_local_manifest_path <- function() {
  env_path <- Sys.getenv("PHARMACOGX_LOCAL_PSET_MANIFEST", "")
  if (nzchar(env_path)) {
    return(path.expand(env_path))
  }

  mcp_dir <- Sys.getenv("PHARMACOGX_MCP_DIR", "")
  candidates <- c(
    if (nzchar(mcp_dir)) {
      file.path(mcp_dir, "local-pset-manifest.csv")
    },
    file.path(getwd(), "inst", "mcp", "local-pset-manifest.csv"),
    file.path(getwd(), "local-pset-manifest.csv")
  )
  candidates <- candidates[nzchar(candidates)]
  existing <- candidates[file.exists(candidates)]
  if (length(existing) > 0) {
    return(existing[[1]])
  }

  candidates[[1]]
}

pgx_default_local_pset_roots <- function() {
  configured <- Sys.getenv("PHARMACOGX_LOCAL_PSET_DIRS", "")
  configured <- if (nzchar(configured)) {
    strsplit(configured, .Platform$path.sep, fixed = TRUE)[[1]]
  } else {
    character()
  }

  package_source <- Sys.getenv("PHARMACOGX_MCP_PACKAGE_SOURCE", "")
  candidates <- c(
    configured,
    if (nzchar(package_source) && !identical(package_source, "installed")) {
      file.path(package_source, "local_data")
    },
    file.path(getwd(), "local_data"),
    file.path(tools::R_user_dir("PharmacoGx", which = "data"), "psets")
  )
  candidates <- unique(normalizePath(
    path.expand(candidates[nzchar(candidates)]),
    mustWork = FALSE
  ))
  candidates[dir.exists(candidates)]
}

pgx_local_alias_from_path <- function(path) {
  stem <- tools::file_path_sans_ext(basename(path))
  known <- c(
    "NCI_ALMANAC_2017" = "nci_almanac",
    "PSet_PRISM" = "prism"
  )
  if (stem %in% names(known)) {
    return(unname(known[[stem]]))
  }

  stem <- sub("^PSet[_-]?", "", stem, ignore.case = TRUE)
  stem <- gsub("[^A-Za-z0-9]+", "_", stem)
  tolower(gsub("(^_+|_+$)", "", stem))
}

pgx_empty_local_manifest <- function() {
  data.frame(
    name = character(),
    path = character(),
    object_type = character(),
    analysis_mode = character(),
    recommended_demo_use = character(),
    aliases = character(),
    notes = character(),
    path_expanded = character(),
    path_normalized = character(),
    exists = logical(),
    extension = character(),
    supported = logical(),
    source = character(),
    size_bytes = numeric(),
    modified_at = character(),
    stringsAsFactors = FALSE
  )
}

pgx_is_likely_pset_path <- function(path) {
  filename <- basename(path)
  directories <- strsplit(
    normalizePath(dirname(path), mustWork = FALSE),
    .Platform$file.sep,
    fixed = TRUE
  )[[1]]
  grepl("pset", filename, ignore.case = TRUE) ||
    tolower(filename) %in% c("gdsc2.rds", "nci_almanac_2017.rds") ||
    any(tolower(directories) == "psets")
}

pgx_discover_local_psets <- function(roots = NULL, recursive = TRUE) {
  if (pgx_is_blank(roots)) {
    roots <- pgx_default_local_pset_roots()
  }
  roots <- pgx_clean_vector(roots)
  if (length(roots) == 0) {
    return(pgx_empty_local_manifest())
  }
  roots <- unique(normalizePath(
    path.expand(roots),
    mustWork = FALSE
  ))
  roots <- roots[dir.exists(roots)]
  if (length(roots) == 0) {
    return(pgx_empty_local_manifest())
  }

  paths <- unique(unlist(
    lapply(
      roots,
      list.files,
      pattern = "\\.(rds|qs)$",
      full.names = TRUE,
      recursive = isTRUE(recursive),
      ignore.case = TRUE
    ),
    use.names = FALSE
  ))
  paths <- paths[vapply(paths, pgx_is_likely_pset_path, logical(1))]
  if (length(paths) == 0) {
    return(pgx_empty_local_manifest())
  }

  info <- file.info(paths)
  data.frame(
    name = vapply(paths, pgx_local_alias_from_path, character(1)),
    path = paths,
    object_type = "PharmacoSet candidate",
    analysis_mode = "local PSet-backed analysis",
    recommended_demo_use = "inspect requirements before analysis",
    aliases = "",
    notes = "Auto-discovered from a configured local PSet root.",
    path_expanded = paths,
    path_normalized = normalizePath(paths, mustWork = FALSE),
    exists = TRUE,
    extension = tolower(tools::file_ext(paths)),
    supported = TRUE,
    source = "discovered",
    size_bytes = unname(info$size),
    modified_at = format(info$mtime, "%Y-%m-%dT%H:%M:%S%z"),
    stringsAsFactors = FALSE
  )
}

pgx_read_local_pset_manifest <- function(
  manifest_path = NULL,
  include_missing = TRUE
) {
  if (pgx_is_blank(manifest_path)) {
    manifest_path <- pgx_default_local_manifest_path()
  }
  manifest_path <- path.expand(manifest_path)

  if (!file.exists(manifest_path)) {
    return(pgx_empty_local_manifest())
  }

  manifest <- utils::read.csv(
    manifest_path,
    stringsAsFactors = FALSE,
    na.strings = c("", "NA")
  )
  required <- c("name", "path")
  missing <- setdiff(required, colnames(manifest))
  if (length(missing) > 0) {
    stop(
      "Local PSet manifest is missing required column(s): ",
      paste(missing, collapse = ", ")
    )
  }

  optional <- setdiff(colnames(pgx_empty_local_manifest()), colnames(manifest))
  optional <- setdiff(
    optional,
    c(
      "path_expanded",
      "path_normalized",
      "exists",
      "extension",
      "supported",
      "source",
      "size_bytes",
      "modified_at"
    )
  )
  for (field in optional) {
    manifest[[field]] <- NA_character_
  }

  relative <- !grepl("^(/|[A-Za-z]:[/\\\\])", manifest$path)
  path_values <- manifest$path
  path_values[relative] <- file.path(
    dirname(normalizePath(manifest_path, mustWork = FALSE)),
    path_values[relative]
  )
  manifest$path_expanded <- path.expand(path_values)
  manifest$exists <- file.exists(manifest$path_expanded)
  manifest$path_normalized <- ifelse(
    manifest$exists,
    normalizePath(manifest$path_expanded, mustWork = FALSE),
    manifest$path_expanded
  )
  manifest$extension <- tolower(tools::file_ext(manifest$path_expanded))
  manifest$supported <- manifest$extension %in% c("qs", "rds")
  manifest$source <- "manifest"
  info <- file.info(manifest$path_expanded)
  manifest$size_bytes <- unname(info$size)
  manifest$modified_at <- format(info$mtime, "%Y-%m-%dT%H:%M:%S%z")

  if (!isTRUE(include_missing)) {
    manifest <- manifest[manifest$exists & manifest$supported, , drop = FALSE]
  }

  manifest
}

pgx_local_pset_catalog <- function(
  manifest_path = NULL,
  roots = NULL,
  discover = TRUE,
  include_missing = TRUE
) {
  manifest <- pgx_read_local_pset_manifest(
    manifest_path = manifest_path,
    include_missing = include_missing
  )
  discovered <- if (isTRUE(discover)) {
    pgx_discover_local_psets(roots = roots)
  } else {
    pgx_empty_local_manifest()
  }
  catalog <- rbind(manifest, discovered)
  if (nrow(catalog) == 0) {
    return(catalog)
  }

  catalog <- catalog[
    !duplicated(catalog$path_normalized),
    ,
    drop = FALSE
  ]
  duplicate_names <- duplicated(catalog$name) |
    duplicated(
      catalog$name,
      fromLast = TRUE
    )
  if (any(duplicate_names)) {
    catalog$name[duplicate_names] <- make.unique(catalog$name[duplicate_names])
  }
  rownames(catalog) <- NULL
  catalog
}

pgx_manifest_alias_values <- function(row) {
  aliases <- pgx_clean_vector(strsplit(row$aliases, ";", fixed = TRUE)[[1]])
  unique(c(row$name, aliases))
}

pgx_resolve_local_pset_alias <- function(dataset) {
  if (pgx_is_blank(dataset)) {
    return(NULL)
  }

  manifest <- pgx_local_pset_catalog(include_missing = TRUE)
  if (nrow(manifest) == 0) {
    return(NULL)
  }

  dataset <- as.character(dataset[[1]])
  exact <- manifest$name == dataset
  alias <- vapply(
    seq_len(nrow(manifest)),
    function(i) {
      dataset %in% pgx_manifest_alias_values(manifest[i, , drop = FALSE])
    },
    logical(1)
  )
  matches <- manifest[exact | alias, , drop = FALSE]
  if (nrow(matches) == 0) {
    return(NULL)
  }
  matches <- matches[order(!(matches$name == dataset)), , drop = FALSE]
  resolved <- matches[1, , drop = FALSE]

  if (!isTRUE(resolved$exists)) {
    stop(
      "Local PSet alias '",
      dataset,
      "' resolves to a missing file: ",
      resolved$path_expanded
    )
  }
  if (!isTRUE(resolved$supported)) {
    stop(
      "Local PSet alias '",
      dataset,
      "' resolves to an unsupported file extension: .",
      resolved$extension
    )
  }

  resolved
}

pgx_is_local_dataset <- function(dataset) {
  if (pgx_is_blank(dataset)) {
    return(FALSE)
  }

  file.exists(path.expand(dataset)) && !dir.exists(path.expand(dataset))
}

pgx_read_local_pset <- function(path) {
  path <- normalizePath(path.expand(path), mustWork = TRUE)
  extension <- tolower(tools::file_ext(path))

  pset <- switch(
    extension,
    "qs" = {
      if (requireNamespace("qs2", quietly = TRUE)) {
        qs2::qs_read(path)
      } else if (requireNamespace("qs", quietly = TRUE)) {
        qs::qread(path)
      } else {
        stop("Reading .qs PharmacoSets requires the qs2 or qs package.")
      }
    },
    "rds" = readRDS(path),
    stop("Unsupported local PSet file extension: .", extension)
  )

  if (!inherits(pset, "PharmacoSet")) {
    stop(
      "Local file did not contain a PharmacoGx PharmacoSet object: ",
      path
    )
  }

  pset
}

pgx_cached_dataset <- function(cache_key, fingerprint, loader) {
  cached <- pgx_dataset_cache[[cache_key]]
  if (!is.null(cached) && identical(cached$fingerprint, fingerprint)) {
    return(cached$pset)
  }

  pset <- loader()
  pgx_dataset_cache[[cache_key]] <- list(
    fingerprint = fingerprint,
    loaded_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    pset = pset
  )
  pset
}

pgx_load_local_pset_cached <- function(path) {
  path <- normalizePath(path.expand(path), mustWork = TRUE)
  info <- file.info(path)
  fingerprint <- paste(
    as.numeric(info$size),
    as.numeric(info$mtime),
    sep = ":"
  )
  pgx_cached_dataset(
    cache_key = paste0("local:", path),
    fingerprint = fingerprint,
    loader = function() pgx_read_local_pset(path)
  )
}

pgx_load_dataset <- function(dataset) {
  if (inherits(dataset, "PharmacoSet")) {
    stop(
      "dataset must be a dataset name, local alias, or .rds/.qs path; ",
      "do not pass a loaded PharmacoSet object to an MCP wrapper."
    )
  }
  if (!is.character(dataset) || length(dataset) != 1 || pgx_is_blank(dataset)) {
    stop("dataset must be one non-empty dataset name, alias, or local path.")
  }

  if (dataset %in% pgx_supported_datasets) {
    return(pgx_cached_dataset(
      cache_key = paste0("bundled:", dataset),
      fingerprint = as.character(utils::packageVersion("PharmacoGx")),
      loader = function() {
        env <- new.env(parent = emptyenv())
        utils::data(list = dataset, package = "PharmacoGx", envir = env)
        env[[dataset]]
      }
    ))
  }

  resolved <- pgx_resolve_local_pset_alias(dataset)
  if (!is.null(resolved)) {
    return(pgx_load_local_pset_cached(resolved$path_normalized))
  }

  if (pgx_is_local_dataset(dataset)) {
    return(pgx_load_local_pset_cached(dataset))
  }

  if (!dataset %in% pgx_supported_datasets) {
    stop(
      "Unsupported dataset '",
      dataset,
      "'. Use one of ",
      paste(pgx_supported_datasets, collapse = ", "),
      ", a local manifest alias, or an explicit local .qs/.rds PharmacoSet file path."
    )
  }
}

pgx_local_dataset_path <- function(dataset) {
  if (!is.character(dataset) || length(dataset) != 1 || pgx_is_blank(dataset)) {
    return(NULL)
  }
  if (pgx_is_local_dataset(dataset)) {
    return(normalizePath(path.expand(dataset), mustWork = TRUE))
  }

  resolved <- tryCatch(
    pgx_resolve_local_pset_alias(dataset),
    error = function(e) NULL
  )
  if (is.null(resolved)) {
    return(NULL)
  }
  normalizePath(resolved$path_normalized, mustWork = TRUE)
}

pgx_dataset_fingerprint <- function(path) {
  path <- normalizePath(path.expand(path), mustWork = TRUE)
  info <- file.info(path)
  paste(
    normalizePath(path, mustWork = TRUE),
    as.numeric(info$size),
    as.numeric(info$mtime),
    sep = "::"
  )
}

pgx_current_tool_file <- function() {
  configured <- Sys.getenv("PHARMACOGX_MCP_TOOL_FILE", "")
  candidates <- c(
    configured,
    file.path(Sys.getenv("PHARMACOGX_MCP_DIR", ""), "PharmacoGx-tools.R"),
    system.file("mcp", "PharmacoGx-tools.R", package = "PharmacoGx"),
    file.path(getwd(), "inst", "mcp", "PharmacoGx-tools.R")
  )
  candidates <- unique(candidates[nzchar(candidates)])
  existing <- candidates[file.exists(candidates)]
  if (length(existing) == 0) {
    stop("Could not locate PharmacoGx-tools.R for isolated capability probing.")
  }
  normalizePath(existing[[1]], mustWork = TRUE)
}

pgx_isolated_probe_script <- function() {
  candidates <- c(
    file.path(Sys.getenv("PHARMACOGX_MCP_DIR", ""), "probe-local-pset.R"),
    system.file("mcp", "probe-local-pset.R", package = "PharmacoGx"),
    file.path(getwd(), "inst", "mcp", "probe-local-pset.R")
  )
  candidates <- unique(candidates[nzchar(candidates)])
  existing <- candidates[file.exists(candidates)]
  if (length(existing) == 0) {
    stop("Could not locate probe-local-pset.R for isolated capability probing.")
  }
  normalizePath(existing[[1]], mustWork = TRUE)
}

pgx_capability_cache_dir <- function() {
  configured <- Sys.getenv("PHARMACOGX_MCP_CAPABILITY_CACHE_DIR", "")
  path <- if (nzchar(configured)) {
    path.expand(configured)
  } else {
    file.path(
      tools::R_user_dir("PharmacoGx", which = "cache"),
      "mcp-capabilities"
    )
  }
  normalizePath(path, mustWork = FALSE)
}

pgx_cache_hash <- function(x) {
  values <- as.integer(charToRaw(enc2utf8(paste(x, collapse = "::"))))
  if (length(values) == 0) {
    return("00000000")
  }
  weights <- (seq_along(values) %% 104729) + 1
  value <- sum((values + 1) * weights) %% 2147483647
  sprintf("%08x", as.integer(value))
}

pgx_capability_cache_file <- function(path, fingerprint, tool_file) {
  stem <- pgx_safe_filename(tools::file_path_sans_ext(basename(path)))
  tool_modified <- as.numeric(file.info(tool_file)$mtime)
  file.path(
    pgx_capability_cache_dir(),
    paste0(stem, "-", pgx_cache_hash(c(fingerprint, tool_modified)), ".rds")
  )
}

pgx_compact_annotation_info <- function(x, patterns) {
  x <- tryCatch(pgx_as_data_frame(x), error = function(e) data.frame())
  if (nrow(x) == 0) {
    return(data.frame())
  }
  fields <- grep(
    paste(patterns, collapse = "|"),
    colnames(x),
    ignore.case = TRUE,
    value = TRUE
  )
  compact <- if (length(fields) == 0) {
    data.frame(row.names = rownames(x))
  } else {
    x[, fields, drop = FALSE]
  }
  for (field in colnames(compact)) {
    if (is.list(compact[[field]])) {
      compact[[field]] <- vapply(
        compact[[field]],
        function(value) paste(pgx_clean_vector(value), collapse = ";"),
        character(1)
      )
    } else {
      compact[[field]] <- pgx_normalize_text(compact[[field]])
    }
  }
  compact
}

pgx_build_dataset_capability_snapshot <- function(dataset) {
  pset <- pgx_load_dataset(dataset)
  treatment_info <- tryCatch(pgx_treatment_info(pset), error = function(e) {
    data.frame()
  })
  if (nrow(treatment_info) == 0) {
    treatment_names <- tryCatch(pgx_treatment_names(pset), error = function(e) {
      character()
    })
    treatment_info <- data.frame(
      treatmentid = treatment_names,
      row.names = treatment_names,
      stringsAsFactors = FALSE
    )
  }
  sample_info <- tryCatch(pgx_sample_info(pset), error = function(e) {
    data.frame()
  })
  if (nrow(sample_info) == 0) {
    sample_names <- tryCatch(pgx_sample_names(pset), error = function(e) {
      character()
    })
    sample_info <- data.frame(
      sampleid = sample_names,
      row.names = sample_names,
      stringsAsFactors = FALSE
    )
  }
  capabilities <- tryCatch(
    pgx_molecular_profile_capabilities(pset),
    error = function(e) data.frame()
  )
  usable_profiles <- if (nrow(capabilities) == 0) {
    character()
  } else {
    capabilities$molecular_profile[capabilities$usable]
  }
  feature_info <- stats::setNames(
    lapply(usable_profiles, function(profile) {
      info <- tryCatch(pgx_feature_info(pset, profile), error = function(e) {
        data.frame()
      })
      pgx_compact_annotation_info(
        info,
        c(
          "feature",
          "gene",
          "symbol",
          "ensembl",
          "probe",
          "id",
          "name",
          "alias"
        )
      )
    }),
    usable_profiles
  )

  local_path <- pgx_local_dataset_path(dataset)
  list(
    schema_version = 1L,
    generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    dataset_path = if (is.null(local_path)) NA_character_ else local_path,
    fingerprint = if (is.null(local_path)) {
      NA_character_
    } else {
      pgx_dataset_fingerprint(local_path)
    },
    object_class = class(pset),
    treatment_info = pgx_compact_annotation_info(
      treatment_info,
      c(
        "treat",
        "drug",
        "compound",
        "name",
        "id",
        "pubchem",
        "cid",
        "inchikey",
        "chembl",
        "alias",
        "brand",
        "target",
        "mechanism"
      )
    ),
    sample_info = pgx_compact_annotation_info(
      sample_info,
      c(
        "sample",
        "cell",
        "name",
        "id",
        "tissue",
        "site",
        "histology",
        "lineage",
        "disease",
        "cancer",
        "origin",
        "panel"
      )
    ),
    sensitivity_measures = tryCatch(
      pgx_sensitivity_measures(pset),
      error = function(e) character()
    ),
    molecular_capabilities = capabilities,
    feature_info = feature_info
  )
}

pgx_tail_file <- function(path, lines = 20L) {
  if (!file.exists(path)) {
    return("")
  }
  values <- readLines(path, warn = FALSE)
  paste(tail(values, lines), collapse = "\n")
}

pgx_probe_local_pset_isolated <- function(path, timeout = NULL) {
  path <- normalizePath(path.expand(path), mustWork = TRUE)
  fingerprint <- pgx_dataset_fingerprint(path)
  tool_file <- pgx_current_tool_file()
  cache_key <- paste(
    fingerprint,
    as.numeric(file.info(tool_file)$mtime),
    sep = "::"
  )
  cached <- pgx_capability_cache[[cache_key]]
  if (!is.null(cached)) {
    return(list(ok = TRUE, snapshot = cached, cache_status = "memory"))
  }

  cache_file <- pgx_capability_cache_file(path, fingerprint, tool_file)
  if (file.exists(cache_file)) {
    cached <- tryCatch(readRDS(cache_file), error = function(e) NULL)
    if (
      is.list(cached) &&
        identical(cached$schema_version, 1L) &&
        identical(cached$fingerprint, fingerprint)
    ) {
      pgx_capability_cache[[cache_key]] <- cached
      return(list(ok = TRUE, snapshot = cached, cache_status = "disk"))
    }
  }

  timeout <- if (pgx_is_blank(timeout)) {
    as.numeric(Sys.getenv("PHARMACOGX_MCP_PROBE_TIMEOUT", "300"))
  } else {
    as.numeric(timeout)
  }
  if (!is.finite(timeout) || timeout <= 0) {
    timeout <- 300
  }

  result_path <- tempfile("pgx-capability-", fileext = ".rds")
  log_path <- tempfile("pgx-capability-", fileext = ".log")
  on.exit(unlink(c(result_path, log_path), force = TRUE), add = TRUE)
  package_source <- Sys.getenv("PHARMACOGX_MCP_PACKAGE_SOURCE", "installed")
  command <- file.path(R.home("bin"), "Rscript")
  status <- tryCatch(
    system2(
      command,
      args = shQuote(c(
        "--vanilla",
        pgx_isolated_probe_script(),
        tool_file,
        package_source,
        path,
        result_path
      )),
      stdout = log_path,
      stderr = log_path,
      timeout = timeout
    ),
    error = function(e) e
  )
  if (inherits(status, "error")) {
    return(list(
      ok = FALSE,
      message = paste("Isolated PSet probe failed:", conditionMessage(status))
    ))
  }
  if (!identical(as.integer(status), 0L) || !file.exists(result_path)) {
    detail <- pgx_tail_file(log_path)
    return(list(
      ok = FALSE,
      message = paste0(
        "Isolated PSet probe exited with status ",
        as.integer(status),
        if (nzchar(detail)) paste0(": ", detail) else ""
      )
    ))
  }

  worker_result <- tryCatch(readRDS(result_path), error = function(e) e)
  if (inherits(worker_result, "error")) {
    return(list(
      ok = FALSE,
      message = paste(
        "Could not read isolated PSet probe result:",
        conditionMessage(worker_result)
      )
    ))
  }
  if (!isTRUE(worker_result$ok)) {
    return(list(ok = FALSE, message = worker_result$message))
  }
  snapshot <- worker_result$snapshot
  if (!identical(snapshot$fingerprint, fingerprint)) {
    return(list(
      ok = FALSE,
      message = "Isolated PSet probe returned stale metadata."
    ))
  }

  dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE)
  cache_temp <- tempfile(
    "pgx-cache-",
    tmpdir = dirname(cache_file),
    fileext = ".rds"
  )
  saveRDS(snapshot, cache_temp, compress = "gzip")
  if (!file.rename(cache_temp, cache_file)) {
    unlink(cache_temp, force = TRUE)
  }
  pgx_capability_cache[[cache_key]] <- snapshot
  list(ok = TRUE, snapshot = snapshot, cache_status = "generated")
}

pgx_dataset_catalog <- function() {
  bundled <- data.frame(
    dataset = c(
      "GDSCsmall",
      "CCLEsmall",
      "CMAPsmall",
      "local .qs/.rds path",
      "HDAC_genes"
    ),
    object_type = c(
      "PharmacoSet",
      "PharmacoSet",
      "PharmacoSet",
      "PharmacoSet",
      "data.frame"
    ),
    analysis_mode = c(
      "drug sensitivity toy dataset",
      "drug sensitivity toy dataset",
      "drug perturbation toy dataset",
      "user-provided local PharmacoSet",
      "example HDAC inhibitor gene signature"
    ),
    recommended_demo_use = c(
      "summarize sensitivity profiles and rank responders",
      "inspect cross-study sensitivity and molecular profile examples",
      "connectivity and perturbation signature examples",
      "real local PSet-backed biomarker, combo, and cross-PSet analysis",
      "signature lookup with CMAPsmall"
    ),
    stringsAsFactors = FALSE
  )

  local_manifest <- tryCatch(
    pgx_local_pset_catalog(include_missing = TRUE),
    error = function(e) pgx_empty_local_manifest()
  )
  if (nrow(local_manifest) == 0) {
    return(bundled)
  }

  local_entries <- data.frame(
    dataset = local_manifest$name,
    object_type = local_manifest$object_type,
    analysis_mode = local_manifest$analysis_mode,
    recommended_demo_use = local_manifest$recommended_demo_use,
    stringsAsFactors = FALSE
  )
  rbind(bundled, local_entries)
}

pgx_trim <- function(x, limit) {
  if (length(x) <= limit) {
    return(x)
  }

  x[seq_len(limit)]
}

pgx_limit_records <- function(records, limit) {
  if (nrow(records) <= limit) {
    return(records)
  }

  records[seq_len(limit), , drop = FALSE]
}

pgx_clean_vector <- function(x) {
  if (is.null(x) || length(x) == 0) {
    return(NULL)
  }

  x <- as.character(x)
  x[!is.na(x) & nzchar(x)]
}

pgx_is_blank <- function(x) {
  is.null(x) || length(x) == 0 || is.na(x[[1]]) || !nzchar(x[[1]])
}

pgx_matrix_to_records <- function(x, value_name = "value", limit = 50) {
  if (length(dim(x)) != 2) {
    stop("Expected a two-dimensional matrix-like object.")
  }
  if (any(dim(x) == 0)) {
    records <- data.frame(
      treatment = character(),
      sample = character(),
      value = numeric(),
      stringsAsFactors = FALSE
    )
    names(records)[[3]] <- value_name
    return(records)
  }

  records <- as.data.frame(as.table(x), stringsAsFactors = FALSE)
  names(records) <- c("treatment", "sample", value_name)
  records <- records[!is.na(records[[value_name]]), , drop = FALSE]
  pgx_limit_records(records, limit)
}

pgx_table_with_id <- function(x, id_name) {
  x <- as.data.frame(x)
  data.frame(
    stats::setNames(list(rownames(x)), id_name),
    x,
    row.names = NULL,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

pgx_plain_matrix <- function(x) {
  if (inherits(x, "SummarizedExperiment")) {
    x <- SummarizedExperiment::assay(x)
  }

  as.matrix(x)
}

pgx_legacy_attr <- function(pset, names) {
  for (name in names) {
    value <- attr(pset, name, exact = TRUE)
    if (!is.null(value)) {
      return(value)
    }
  }

  NULL
}

pgx_has_legacy_attrs <- function(pset) {
  any(
    !vapply(
      c("cell", "drug", "sensitivity"),
      function(name) is.null(attr(pset, name, exact = TRUE)),
      logical(1)
    )
  )
}

pgx_as_data_frame <- function(x) {
  if (is.null(x)) {
    return(data.frame())
  }

  as.data.frame(x, stringsAsFactors = FALSE)
}

pgx_sample_info <- function(pset) {
  if (pgx_has_legacy_attrs(pset)) {
    return(pgx_as_data_frame(pgx_legacy_attr(pset, c("sample", "cell"))))
  }

  info <- tryCatch(
    suppressWarnings(PharmacoGx::sampleInfo(pset)),
    error = function(e) NULL
  )
  if (!is.null(info)) {
    return(pgx_as_data_frame(info))
  }

  pgx_as_data_frame(pgx_legacy_attr(pset, c("sample", "cell")))
}

pgx_treatment_info <- function(pset) {
  if (pgx_has_legacy_attrs(pset)) {
    return(pgx_as_data_frame(pgx_legacy_attr(pset, c("treatment", "drug"))))
  }

  info <- tryCatch(
    suppressWarnings(PharmacoGx::treatmentInfo(pset)),
    error = function(e) NULL
  )
  if (!is.null(info)) {
    return(pgx_as_data_frame(info))
  }

  pgx_as_data_frame(pgx_legacy_attr(pset, c("treatment", "drug")))
}

pgx_id_values <- function(info, preferred_columns) {
  if (nrow(info) == 0) {
    return(character())
  }

  values <- character()
  for (field in intersect(preferred_columns, colnames(info))) {
    values <- c(values, as.character(info[[field]]))
  }
  row_ids <- rownames(info)
  if (
    !is.null(row_ids) && !identical(row_ids, as.character(seq_len(nrow(info))))
  ) {
    values <- c(row_ids, values)
  }

  unique(pgx_clean_vector(values))
}

pgx_sample_names <- function(pset) {
  if (pgx_has_legacy_attrs(pset)) {
    return(pgx_id_values(
      pgx_sample_info(pset),
      c("sampleid", "cellid", "sample", "cell", "depmap_id", "COSMICID")
    ))
  }

  values <- tryCatch(
    suppressWarnings(PharmacoGx::sampleNames(pset)),
    error = function(e) NULL
  )
  if (!is.null(values) && length(values) > 0) {
    return(as.character(values))
  }

  pgx_id_values(
    pgx_sample_info(pset),
    c("sampleid", "cellid", "sample", "cell", "depmap_id", "COSMICID")
  )
}

pgx_treatment_names <- function(pset) {
  if (pgx_has_legacy_attrs(pset)) {
    return(pgx_id_values(
      pgx_treatment_info(pset),
      c("treatmentid", "drugid", "treatment", "drug", "name", "compound")
    ))
  }

  values <- tryCatch(
    suppressWarnings(PharmacoGx::treatmentNames(pset)),
    error = function(e) NULL
  )
  if (!is.null(values) && length(values) > 0) {
    return(as.character(values))
  }

  pgx_id_values(
    pgx_treatment_info(pset),
    c("treatmentid", "drugid", "treatment", "drug", "name", "compound")
  )
}

pgx_molecular_profiles <- function(pset) {
  if (pgx_has_legacy_attrs(pset)) {
    profiles <- pgx_legacy_attr(pset, "molecularProfiles")
    if (is.null(profiles)) {
      return(list())
    }
    return(profiles)
  }

  profiles <- tryCatch(
    suppressWarnings(PharmacoGx::molecularProfilesSlot(pset)),
    error = function(e) NULL
  )
  if (!is.null(profiles)) {
    return(profiles)
  }

  profiles <- pgx_legacy_attr(pset, "molecularProfiles")
  if (is.null(profiles)) {
    list()
  } else {
    profiles
  }
}

pgx_mdata_names <- function(pset) {
  values <- tryCatch(
    suppressWarnings(PharmacoGx::mDataNames(pset)),
    error = function(e) NULL
  )
  if (!is.null(values) && length(values) > 0) {
    return(as.character(values))
  }

  names(pgx_molecular_profiles(pset))
}

pgx_feature_info <- function(pset, mDataType) {
  if (pgx_has_legacy_attrs(pset)) {
    profile <- pgx_molecular_profiles(pset)[[mDataType]]
    if (!is.null(profile)) {
      row_data <- tryCatch(
        as.data.frame(SummarizedExperiment::rowData(profile)),
        error = function(e) data.frame()
      )
      if (nrow(row_data) > 0) {
        return(row_data)
      }

      mat <- pgx_plain_matrix(profile)
      return(data.frame(row.names = rownames(mat)))
    }
  }

  info <- tryCatch(
    suppressWarnings(PharmacoGx::featureInfo(pset, mDataType)),
    error = function(e) NULL
  )
  if (!is.null(info)) {
    return(pgx_as_data_frame(info))
  }

  profile <- pgx_molecular_profiles(pset)[[mDataType]]
  if (is.null(profile)) {
    return(data.frame())
  }

  row_data <- tryCatch(
    as.data.frame(SummarizedExperiment::rowData(profile)),
    error = function(e) data.frame()
  )
  if (nrow(row_data) > 0) {
    return(row_data)
  }

  mat <- pgx_plain_matrix(profile)
  data.frame(row.names = rownames(mat))
}

pgx_molecular_profile_capabilities <- function(pset) {
  profile_names <- pgx_mdata_names(pset)
  if (length(profile_names) == 0) {
    return(data.frame(
      molecular_profile = character(),
      usable = logical(),
      feature_count = integer(),
      sample_count = integer(),
      assay_names = character(),
      feature_fields = character(),
      stringsAsFactors = FALSE
    ))
  }

  profiles <- pgx_molecular_profiles(pset)
  rows <- lapply(profile_names, function(profile_name) {
    profile <- profiles[[profile_name]]
    dims <- tryCatch(dim(profile), error = function(e) NULL)
    if (is.null(dims) || length(dims) < 2) {
      mat <- tryCatch(pgx_plain_matrix(profile), error = function(e) NULL)
      dims <- if (is.null(mat)) c(0L, 0L) else dim(mat)
    }
    feature_info <- tryCatch(
      pgx_feature_info(pset, profile_name),
      error = function(e) data.frame()
    )
    assay_names <- tryCatch(
      SummarizedExperiment::assayNames(profile),
      error = function(e) character()
    )
    feature_count <- if (length(dims) >= 1) as.integer(dims[[1]]) else 0L
    sample_count <- if (length(dims) >= 2) as.integer(dims[[2]]) else 0L
    data.frame(
      molecular_profile = profile_name,
      usable = feature_count > 0L && sample_count > 0L,
      feature_count = feature_count,
      sample_count = sample_count,
      assay_names = paste(assay_names, collapse = ";"),
      feature_fields = paste(colnames(feature_info), collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

pgx_summarize_profile_by_sample <- function(
  profile,
  features,
  samples,
  summary_stat
) {
  summary_stat <- match.arg(
    summary_stat,
    c("mean", "median", "first", "last", "and", "or")
  )
  feature_ids <- intersect(features, rownames(profile))
  sample_ids <- colnames(profile)
  if (inherits(profile, "SummarizedExperiment")) {
    annotation <- as.data.frame(SummarizedExperiment::colData(profile))
    id_fields <- intersect(c("cellid", "sampleid"), names(annotation))
    if (length(id_fields) > 0L) {
      sample_ids <- as.character(annotation[[id_fields[[1L]]]])
    }
  }
  if (
    length(sample_ids) != ncol(profile) ||
      anyNA(sample_ids) ||
      any(!nzchar(trimws(sample_ids)))
  ) {
    stop("Molecular assay columns have missing biological sample identifiers.")
  }
  sample_keep <- intersect(samples, sample_ids)
  # Subset features before materializing a potentially large expression assay.
  mat <- pgx_plain_matrix(profile[feature_ids, , drop = FALSE])
  result <- matrix(
    NA_real_,
    nrow = length(feature_ids),
    ncol = length(sample_keep),
    dimnames = list(feature_ids, sample_keep)
  )
  summarize_values <- function(values) {
    if (summary_stat == "first") {
      return(values[[1L]])
    }
    if (summary_stat == "last") {
      return(values[[length(values)]])
    }
    if (summary_stat %in% c("and", "or")) {
      if (!all(is.na(values) | values %in% c(0, 1))) {
        stop("Logical molecular aggregation requires binary 0/1 values.")
      }
      if (all(is.na(values))) {
        return(NA_real_)
      }
      return(as.numeric(
        if (summary_stat == "and") {
          all(as.logical(as.numeric(values)), na.rm = TRUE)
        } else {
          any(as.logical(as.numeric(values)), na.rm = TRUE)
        }
      ))
    }
    if (!is.numeric(values)) {
      stop("Mean/median molecular aggregation requires numeric values.")
    }
    values <- values[is.finite(values)]
    if (length(values) == 0L) {
      return(NA_real_)
    }
    if (summary_stat == "median") median(values) else mean(values)
  }
  for (sample in sample_keep) {
    columns <- which(sample_ids == sample)
    result[, sample] <- vapply(
      seq_along(feature_ids),
      function(i) {
        summarize_values(mat[i, columns])
      },
      numeric(1)
    )
  }
  result
}

pgx_summarize_molecular_profiles <- function(
  pset,
  mDataType,
  features,
  cell.lines,
  summary.stat
) {
  if (pgx_has_legacy_attrs(pset)) {
    profile <- pgx_molecular_profiles(pset)[[mDataType]]
    if (is.null(profile)) {
      stop("Molecular profile '", mDataType, "' was not found.")
    }

    return(pgx_summarize_profile_by_sample(
      profile,
      features,
      cell.lines,
      summary.stat
    ))
  }

  summarized <- tryCatch(
    suppressWarnings(PharmacoGx::summarizeMolecularProfiles(
      pset,
      mDataType = mDataType,
      features = features,
      cell.lines = cell.lines,
      summary.stat = summary.stat,
      verbose = FALSE
    )),
    error = function(e) NULL
  )
  if (!is.null(summarized)) {
    return(pgx_plain_matrix(summarized))
  }

  profile <- pgx_molecular_profiles(pset)[[mDataType]]
  if (is.null(profile)) {
    stop("Molecular profile '", mDataType, "' was not found.")
  }

  pgx_summarize_profile_by_sample(profile, features, cell.lines, summary.stat)
}

pgx_select_fields <- function(x, id_name, fields = NULL) {
  x <- pgx_table_with_id(x, id_name)
  fields <- pgx_clean_vector(fields)

  if (is.null(fields)) {
    return(x)
  }

  missing_fields <- setdiff(fields, colnames(x))
  if (length(missing_fields) > 0) {
    stop(
      "Unknown field(s): ",
      paste(missing_fields, collapse = ", "),
      ". Available fields are: ",
      paste(colnames(x), collapse = ", ")
    )
  }

  keep <- unique(c(id_name, fields))
  x[, keep, drop = FALSE]
}

pgx_resolve_match_mode <- function(match_mode) {
  match.arg(match_mode, c("exact", "contains"))
}

pgx_candidate_columns <- function(x, patterns) {
  candidates <- unique(c(
    rownames(x),
    grep(
      paste(patterns, collapse = "|"),
      colnames(x),
      ignore.case = TRUE,
      value = TRUE
    )
  ))
  intersect(candidates, colnames(x))
}

pgx_find_rows <- function(x, queries, match_mode, search_fields) {
  queries <- pgx_clean_vector(queries)
  if (is.null(queries)) {
    return(integer())
  }

  match_mode <- pgx_resolve_match_mode(match_mode)
  searchable <- pgx_table_with_id(x, ".row_id")
  search_fields <- unique(c(".row_id", search_fields))
  search_fields <- intersect(search_fields, colnames(searchable))

  matched <- logical(nrow(searchable))
  for (query in queries) {
    if (identical(match_mode, "exact")) {
      field_match <- Reduce(
        `|`,
        lapply(search_fields, function(field) {
          pgx_casefold(searchable[[field]]) == pgx_casefold(query)
        })
      )
    } else {
      field_match <- Reduce(
        `|`,
        lapply(search_fields, function(field) {
          grepl(
            pgx_casefold(query),
            pgx_casefold(searchable[[field]]),
            fixed = TRUE
          )
        })
      )
    }
    matched <- matched | field_match
  }

  which(matched)
}

pgx_matched_fields <- function(x, rows, queries, match_mode, search_fields) {
  if (length(rows) == 0) {
    return(character())
  }
  searchable <- pgx_table_with_id(x, ".row_id")
  search_fields <- intersect(
    unique(c(".row_id", search_fields)),
    colnames(searchable)
  )
  queries <- pgx_clean_vector(queries)
  match_mode <- pgx_resolve_match_mode(match_mode)

  vapply(
    rows,
    function(row_index) {
      matched <- vapply(
        search_fields,
        function(field) {
          values <- pgx_casefold(searchable[[field]][[row_index]])
          any(vapply(
            queries,
            function(query) {
              query <- pgx_casefold(query)
              if (identical(match_mode, "exact")) {
                values == query
              } else {
                grepl(query, values, fixed = TRUE)
              }
            },
            logical(1)
          ))
        },
        logical(1)
      )
      paste(search_fields[matched], collapse = ";")
    },
    character(1)
  )
}

pgx_normalize_identifier_values <- function(x) {
  x <- pgx_normalize_text(x)
  if (length(x) == 0) {
    return(character())
  }
  x <- trimws(x)
  x <- sub(
    paste0(
      "^(pubchem([ _-]*(cid|compound))?|cid|chembl([ _-]*id)?|",
      "ensembl([ _-]*(gene|transcript|protein))?([ _-]*id)?|",
      "entrez([ _-]*(gene))?([ _-]*id)?)\\s*[:=#]?\\s*"
    ),
    "",
    x,
    ignore.case = TRUE,
    perl = TRUE
  )
  x <- sub(
    "^(ENS[A-Z]*[GTP][0-9]+)\\.[0-9]+$",
    "\\1",
    x,
    ignore.case = TRUE,
    perl = TRUE
  )
  sub("^([0-9]+)\\.0+$", "\\1", x)
}

pgx_identifier_tokens <- function(x) {
  x <- pgx_clean_vector(x)
  if (is.null(x)) {
    return(character())
  }

  tokens <- unlist(
    strsplit(x, "\\s*(?:///|;|\\|)\\s*", perl = TRUE),
    use.names = FALSE
  )
  unique(pgx_clean_vector(pgx_normalize_identifier_values(tokens)))
}

pgx_identifier_matches <- function(values, query_tokens) {
  values <- pgx_normalize_text(values)
  normalized <- pgx_casefold(pgx_normalize_identifier_values(values))
  matched <- normalized %in% query_tokens
  delimited <- grepl("///|;|\\|", values)
  if (any(delimited)) {
    matched[delimited] <- vapply(
      values[delimited],
      function(value) {
        value_tokens <- pgx_casefold(pgx_identifier_tokens(value))
        length(intersect(query_tokens, value_tokens)) > 0
      },
      logical(1)
    )
  }
  matched
}

pgx_resolve_annotation_queries <- function(
  x,
  queries,
  search_fields = colnames(x)
) {
  queries <- pgx_clean_vector(queries)
  if (is.null(queries) || nrow(x) == 0) {
    return(list(
      rows = integer(),
      matched_queries = character(),
      unmatched_queries = if (is.null(queries)) character() else queries,
      matched_fields = character()
    ))
  }

  searchable <- pgx_table_with_id(x, ".row_id")
  search_fields <- intersect(
    unique(c(".row_id", search_fields)),
    colnames(searchable)
  )
  query_rows <- stats::setNames(vector("list", length(queries)), queries)
  query_fields <- stats::setNames(vector("list", length(queries)), queries)

  for (query in queries) {
    query_tokens <- pgx_casefold(pgx_identifier_tokens(query))
    for (field in search_fields) {
      values <- pgx_normalize_text(searchable[[field]])
      matched <- pgx_identifier_matches(values, query_tokens)
      if (any(matched)) {
        query_rows[[query]] <- union(query_rows[[query]], which(matched))
        query_fields[[query]] <- union(query_fields[[query]], field)
      }
    }
  }

  matched <- names(query_rows)[lengths(query_rows) > 0]
  list(
    rows = sort(unique(unlist(query_rows, use.names = FALSE))),
    matched_queries = matched,
    unmatched_queries = setdiff(queries, matched),
    matched_fields = unique(unlist(query_fields[matched], use.names = FALSE))
  )
}

pgx_first_row_identifier <- function(info, rows, preferred_fields) {
  if (length(rows) == 0) {
    return(character())
  }

  searchable <- pgx_table_with_id(info, ".row_id")
  preferred_fields <- intersect(
    c(preferred_fields, ".row_id"),
    colnames(searchable)
  )
  vapply(
    rows,
    function(row) {
      values <- pgx_clean_vector(unlist(
        searchable[row, preferred_fields, drop = FALSE],
        use.names = FALSE
      ))
      if (length(values) == 0) "" else values[[1]]
    },
    character(1)
  )
}

pgx_resolve_treatment_identifiers <- function(pset, queries) {
  info <- pgx_treatment_info(pset)
  if (nrow(info) == 0) {
    treatment_names <- pgx_treatment_names(pset)
    info <- data.frame(
      treatmentid = treatment_names,
      row.names = treatment_names,
      stringsAsFactors = FALSE
    )
  }
  resolution <- pgx_resolve_annotation_queries(
    info,
    queries,
    search_fields = colnames(info)
  )
  resolution$identifiers <- pgx_first_row_identifier(
    info,
    resolution$rows,
    c("treatmentid", "drugid", "treatment", "drug", "name", "compound")
  )
  resolution
}

pgx_sample_context_variants <- function(query) {
  query <- pgx_clean_vector(query)
  if (is.null(query)) {
    return(character())
  }

  aliases <- list(
    pdac = c(
      "PDAC",
      "pancreas",
      "pancreatic",
      "pancreatic adenocarcinoma",
      "pancreatic ductal adenocarcinoma"
    ),
    nsclc = c("NSCLC", "non-small cell lung cancer", "lung"),
    gbm = c("GBM", "glioblastoma", "brain")
  )
  expanded <- unlist(lapply(query, function(value) {
    key <- pgx_casefold(value)
    if (key %in% names(aliases)) aliases[[key]] else value
  }))
  unique(pgx_clean_vector(expanded))
}

pgx_profile_family <- function(x) {
  x <- pgx_casefold(x)
  if (grepl("mirna|micro[ _-]*rna", x)) {
    "mirna"
  } else if (grepl("rna|transcript|expression", x)) {
    "rna"
  } else if (grepl("mutation|variant|snv", x)) {
    "mutation"
  } else if (grepl("cnv|copy[ _-]*number", x)) {
    "cnv"
  } else if (grepl("methyl", x)) {
    "methylation"
  } else if (grepl("protein|proteom", x)) {
    "protein"
  } else {
    x
  }
}

pgx_resolve_molecular_profile <- function(requested, capabilities) {
  if (pgx_is_blank(requested) || nrow(capabilities) == 0) {
    return(capabilities[FALSE, , drop = FALSE])
  }

  exact <- pgx_casefold(capabilities$molecular_profile) ==
    pgx_casefold(requested)
  if (any(exact)) {
    return(capabilities[exact, , drop = FALSE])
  }
  family <- pgx_profile_family(requested)
  candidates <- vapply(
    capabilities$molecular_profile,
    pgx_profile_family,
    character(1)
  ) ==
    family
  capabilities[candidates, , drop = FALSE]
}

pgx_parse_filters <- function(filters) {
  if (is.null(filters) || length(filters) == 0) {
    stop(
      "filters must be a named object or JSON object with at least one field."
    )
  }

  if (is.character(filters) && length(filters) == 1) {
    filters <- jsonlite::fromJSON(filters, simplifyVector = FALSE)
  }

  if (
    !is.list(filters) || is.null(names(filters)) || any(!nzchar(names(filters)))
  ) {
    stop(
      "filters must be a named object, for example {\"tissueid\": \"lung\"}."
    )
  }

  filters
}

pgx_match_filter_values <- function(values, expected, match_mode) {
  expected <- pgx_clean_vector(expected)
  if (is.null(expected)) {
    return(rep(TRUE, length(values)))
  }

  values <- pgx_normalize_text(values)
  match_mode <- pgx_resolve_match_mode(match_mode)

  Reduce(
    `|`,
    lapply(expected, function(value) {
      if (identical(match_mode, "exact")) {
        pgx_casefold(values) == pgx_casefold(value)
      } else {
        grepl(pgx_casefold(value), pgx_casefold(values), fixed = TRUE)
      }
    })
  )
}

pgx_metric_direction <- function(metric, rank_direction = "auto") {
  rank_direction <- match.arg(rank_direction, c("auto", "lowest", "highest"))
  if (!identical(rank_direction, "auto")) {
    return(rank_direction)
  }

  if (grepl("aac|amax", metric, ignore.case = TRUE)) {
    "highest"
  } else {
    "lowest"
  }
}

pgx_metric_suggestions <- function(metric, available) {
  available <- pgx_clean_vector(available)
  if (is.null(available)) {
    return(character())
  }
  metric <- pgx_casefold(metric)
  tokens <- unique(unlist(strsplit(metric, "[_ .-]+")))
  scores <- vapply(
    available,
    function(candidate) {
      candidate_folded <- pgx_casefold(candidate)
      sum(vapply(tokens, grepl, logical(1), x = candidate_folded, fixed = TRUE))
    },
    numeric(1)
  )
  available[order(scores, decreasing = TRUE)][seq_len(min(
    5L,
    length(available)
  ))]
}

pgx_resolve_metric <- function(metric, available) {
  available <- pgx_clean_vector(available)
  if (is.null(available)) {
    return(list(resolved = NULL, suggestions = character()))
  }
  if (pgx_is_blank(metric) || identical(pgx_casefold(metric), "auto")) {
    preferred <- c(
      "aac_recomputed",
      "auc_recomputed",
      "aac_published",
      "auc_published",
      "ic50_recomputed",
      "ic50_published"
    )
    resolved <- intersect(preferred, available)
    if (length(resolved) == 0) {
      resolved <- available
    }
    return(list(resolved = resolved[[1]], suggestions = available))
  }
  exact <- available[pgx_casefold(available) == pgx_casefold(metric)]
  if (length(exact) > 0) {
    return(list(
      resolved = exact[[1]],
      suggestions = pgx_metric_suggestions(metric, available),
      match_type = "exact"
    ))
  }

  legacy_aliases <- c(
    auc_recomputed = "aac_recomputed",
    auc_recomputed_star = "aac_recomputed_star"
  )
  metric_folded <- pgx_casefold(metric)
  alias <- unname(legacy_aliases[metric_folded])
  alias <- alias[!is.na(alias)]
  if (length(alias) == 0) {
    reverse_alias <- names(legacy_aliases)[
      pgx_casefold(legacy_aliases) == metric_folded
    ]
    if (length(reverse_alias) > 0) {
      alias <- reverse_alias[[1]]
    }
  }
  alias_match <- if (length(alias) == 0) {
    character()
  } else {
    available[pgx_casefold(available) == pgx_casefold(alias)]
  }
  list(
    resolved = if (length(alias_match) > 0) alias_match[[1]] else NULL,
    suggestions = pgx_metric_suggestions(metric, available),
    match_type = if (length(alias_match) > 0) "legacy_alias" else "unresolved"
  )
}

pgx_first_present_column <- function(x, candidates) {
  matched <- intersect(candidates, colnames(x))
  if (length(matched) == 0) {
    return(NULL)
  }

  matched[[1]]
}

pgx_sensitivity_id_columns <- c(
  "treatmentid",
  "treatment1id",
  "drugid",
  "drug",
  "treatment",
  "sampleid",
  "cellid",
  "cell",
  "sample",
  "tech_rep",
  "replicate",
  "dose",
  "treatment1dose",
  "viability",
  "mean_viability"
)

pgx_numeric_metric_columns <- function(x) {
  numeric_cols <- vapply(x, is.numeric, logical(1))
  setdiff(colnames(x)[numeric_cols], pgx_sensitivity_id_columns)
}

pgx_sensitivity_metric_catalog <- function(pset) {
  response <- pgx_treatment_response(pset)
  table_names <- c(
    intersect(c("mono_profiles", "profiles"), names(response)),
    setdiff(names(response), c("mono_profiles", "profiles"))
  )
  rows <- lapply(table_names, function(table_name) {
    x <- response[[table_name]]
    metrics <- pgx_numeric_metric_columns(x)
    if (length(metrics) == 0) {
      return(NULL)
    }
    data.frame(
      table = table_name,
      metric = metrics,
      stringsAsFactors = FALSE
    )
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0) {
    return(data.frame(
      table = character(),
      metric = character(),
      stringsAsFactors = FALSE
    ))
  }

  unique(do.call(rbind, rows))
}

pgx_sensitivity_table_for_metric <- function(pset, sensitivity.measure) {
  response <- pgx_treatment_response(pset)
  preferred <- c(
    intersect(c("mono_profiles", "profiles"), names(response)),
    setdiff(names(response), c("mono_profiles", "profiles"))
  )
  matching <- preferred[vapply(
    preferred,
    function(table_name) {
      sensitivity.measure %in%
        colnames(response[[table_name]]) &&
        is.numeric(response[[table_name]][[sensitivity.measure]])
    },
    logical(1)
  )]
  if (length(matching) == 0) {
    catalog <- pgx_sensitivity_metric_catalog(pset)
    available <- unique(catalog$metric)
    stop(
      "Sensitivity measure '",
      sensitivity.measure,
      "' was not found in a monotherapy profile table. Available measures are: ",
      if (length(available) > 0) paste(available, collapse = ", ") else "none",
      "."
    )
  }

  list(name = matching[[1]], table = response[[matching[[1]]]])
}

pgx_aligned_sensitivity_info <- function(pset, profiles, columns = NULL) {
  info <- pgx_sensitivity_info(pset)
  if (nrow(info) == 0 || nrow(profiles) == 0) {
    return(data.frame())
  }
  if (!is.null(columns)) {
    info <- info[, intersect(columns, colnames(info)), drop = FALSE]
  }

  profile_ids <- rownames(profiles)
  info_ids <- rownames(info)
  if (
    !is.null(profile_ids) &&
      !is.null(info_ids) &&
      !identical(profile_ids, as.character(seq_len(nrow(profiles))))
  ) {
    index <- match(profile_ids, info_ids)
    if (any(!is.na(index))) {
      return(info[index, , drop = FALSE])
    }
  }

  if (nrow(info) == nrow(profiles)) {
    return(info)
  }

  data.frame()
}

pgx_profile_identifier <- function(
  pset,
  profiles,
  profile_candidates,
  info_candidates
) {
  profile_col <- pgx_first_present_column(profiles, profile_candidates)
  if (!is.null(profile_col)) {
    return(list(
      values = as.character(profiles[[profile_col]]),
      source = profile_col
    ))
  }

  info <- pgx_aligned_sensitivity_info(
    pset,
    profiles,
    columns = info_candidates
  )
  info_col <- pgx_first_present_column(info, info_candidates)
  if (!is.null(info_col)) {
    return(list(
      values = as.character(info[[info_col]]),
      source = paste0("sensitivityInfo:", info_col)
    ))
  }

  NULL
}

pgx_resolve_metadata_queries <- function(
  info,
  queries,
  available_values,
  preferred_id_columns
) {
  queries <- pgx_clean_vector(queries)
  available_values <- unique(pgx_clean_vector(available_values))
  if (is.null(queries)) {
    return(data.frame(
      query = character(),
      resolved = character(),
      stringsAsFactors = FALSE
    ))
  }

  rows <- lapply(queries, function(query) {
    direct <- available_values[
      pgx_casefold(available_values) == pgx_casefold(query)
    ]
    if (length(direct) > 0) {
      return(data.frame(
        query = query,
        resolved = direct,
        stringsAsFactors = FALSE
      ))
    }
    if (nrow(info) == 0) {
      return(NULL)
    }

    searchable <- data.frame(
      .row = rownames(info),
      info,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    matches <- apply(
      searchable,
      1,
      function(row) {
        any(pgx_casefold(row) == pgx_casefold(query), na.rm = TRUE)
      }
    )
    if (!any(matches)) {
      return(NULL)
    }

    id_columns <- unique(c(
      intersect(preferred_id_columns, colnames(searchable)),
      colnames(searchable)
    ))
    candidates <- unique(unlist(
      searchable[matches, id_columns, drop = FALSE],
      use.names = FALSE
    ))
    candidates <- pgx_clean_vector(as.character(candidates))
    resolved <- available_values[
      pgx_casefold(available_values) %in% pgx_casefold(candidates)
    ]
    if (length(resolved) == 0) {
      return(NULL)
    }

    data.frame(
      query = query,
      resolved = resolved,
      stringsAsFactors = FALSE
    )
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0) {
    return(data.frame(
      query = character(),
      resolved = character(),
      stringsAsFactors = FALSE
    ))
  }

  unique(do.call(rbind, rows))
}

pgx_resolve_treatment_queries <- function(pset, queries, available_values) {
  pgx_resolve_metadata_queries(
    info = pgx_treatment_info(pset),
    queries = queries,
    available_values = available_values,
    preferred_id_columns = c(
      "treatmentid",
      "drugid",
      "NSC_number",
      "NSC",
      "name",
      "compound",
      ".row"
    )
  )
}

pgx_summary_treatment_rows <- function(pset, summary, query) {
  available <- rownames(summary)
  resolution <- attr(summary, "resolved_treatments", exact = TRUE)
  if (is.null(resolution) || nrow(resolution) == 0) {
    resolution <- pgx_resolve_treatment_queries(pset, query, available)
  }
  unique(resolution$resolved[resolution$query == query])
}

pgx_summarize_sensitivity_profiles <- function(
  pset,
  sensitivity.measure,
  drugs = NULL,
  cell.lines = NULL
) {
  if (pgx_has_legacy_attrs(pset)) {
    return(pgx_summarize_legacy_sensitivity_profiles(
      pset = pset,
      sensitivity.measure = sensitivity.measure,
      drugs = drugs,
      cell.lines = cell.lines
    ))
  }

  args <- list(
    object = pset,
    sensitivity.measure = sensitivity.measure,
    verbose = FALSE
  )
  if (!is.null(drugs)) {
    args$drugs <- drugs
  }
  if (!is.null(cell.lines)) {
    args$cell.lines <- cell.lines
  }

  summary <- tryCatch(
    suppressWarnings(do.call(PharmacoGx::summarizeSensitivityProfiles, args)),
    error = function(e) NULL
  )
  if (!is.null(summary)) {
    return(summary)
  }

  pgx_summarize_legacy_sensitivity_profiles(
    pset = pset,
    sensitivity.measure = sensitivity.measure,
    drugs = drugs,
    cell.lines = cell.lines
  )
}

pgx_summarize_legacy_sensitivity_profiles <- function(
  pset,
  sensitivity.measure,
  drugs = NULL,
  cell.lines = NULL
) {
  selected <- pgx_sensitivity_table_for_metric(pset, sensitivity.measure)
  profiles <- selected$table
  treatment <- pgx_profile_identifier(
    pset = pset,
    profiles = profiles,
    profile_candidates = c(
      "treatmentid",
      "treatment1id",
      "drugid",
      "drug",
      "treatment"
    ),
    info_candidates = c("treatmentid", "drugid", "drug", "treatment")
  )
  sample <- pgx_profile_identifier(
    pset = pset,
    profiles = profiles,
    profile_candidates = c("sampleid", "cellid", "cell", "sample"),
    info_candidates = c("sampleid", "cellid", "cell", "sample")
  )
  if (is.null(treatment) || is.null(sample)) {
    stop(
      "Could not resolve treatment and sample identifiers for response table '",
      selected$name,
      "'."
    )
  }

  values <- profiles[[sensitivity.measure]]
  keep <- !is.na(values) & is.finite(values)
  resolved_treatments <- pgx_resolve_treatment_queries(
    pset,
    drugs,
    unique(treatment$values[keep])
  )
  if (!is.null(drugs)) {
    keep <- keep &
      pgx_casefold(treatment$values) %in%
        pgx_casefold(resolved_treatments$resolved)
  }
  if (!is.null(cell.lines)) {
    keep <- keep & pgx_casefold(sample$values) %in% pgx_casefold(cell.lines)
  }
  x <- data.frame(
    treatment = treatment$values[keep],
    sample = sample$values[keep],
    value = values[keep],
    stringsAsFactors = FALSE
  )
  if (nrow(x) == 0) {
    result <- matrix(numeric(), nrow = 0, ncol = 0)
    attr(result, "response_table") <- selected$name
    attr(result, "resolved_treatments") <- resolved_treatments
    return(result)
  }

  summarized <- stats::aggregate(
    value ~ treatment + sample,
    data = x,
    FUN = median,
    na.rm = TRUE
  )
  treatments <- unique(summarized$treatment)
  samples <- unique(summarized$sample)
  mat <- matrix(
    NA_real_,
    nrow = length(treatments),
    ncol = length(samples),
    dimnames = list(treatments, samples)
  )
  mat[cbind(
    match(summarized$treatment, treatments),
    match(summarized$sample, samples)
  )] <- summarized$value
  attr(mat, "response_table") <- selected$name
  attr(mat, "resolved_treatments") <- resolved_treatments
  attr(mat, "treatment_identifier_source") <- treatment$source
  attr(mat, "sample_identifier_source") <- sample$source
  mat
}

pgx_sensitivity_records <- function(pset, metric, drug = NULL) {
  summary <- pgx_summarize_sensitivity_profiles(
    pset = pset,
    sensitivity.measure = metric,
    drugs = if (!pgx_is_blank(drug)) drug else NULL
  )
  if (length(dim(summary)) != 2 || any(dim(summary) == 0)) {
    records <- data.frame(
      treatment = character(),
      sample = character(),
      value = numeric(),
      stringsAsFactors = FALSE
    )
    names(records)[[3]] <- metric
    return(records)
  }

  records <- as.data.frame(as.table(summary), stringsAsFactors = FALSE)
  names(records) <- c("treatment", "sample", metric)
  records <- records[!is.na(records[[metric]]), , drop = FALSE]
  records
}

pgx_rank_metric_records <- function(records, metric, direction, top_n, suffix) {
  if (nrow(records) == 0) {
    return(data.frame())
  }

  ranked <- records[
    order(
      records[[metric]],
      decreasing = identical(direction, "highest")
    ),
    ,
    drop = FALSE
  ]
  ranked <- ranked[seq_len(min(nrow(ranked), top_n)), , drop = FALSE]
  ranked$rank <- seq_len(nrow(ranked))
  ranked$direction <- direction
  ranked$metric <- metric
  names(ranked)[names(ranked) == metric] <- suffix
  ranked[, c("metric", "rank", "direction", "treatment", "sample", suffix)]
}

pgx_safe_filename <- function(...) {
  x <- paste(..., sep = "_")
  x <- gsub("[^A-Za-z0-9_.-]+", "_", x)
  x <- gsub("_+", "_", x)
  x
}

pgx_resolve_experiments <- function(pset, drug, sample) {
  info <- pgx_sensitivity_info(pset)
  required <- c("sampleid", "treatmentid")
  missing_cols <- setdiff(required, colnames(info))
  if (length(missing_cols) > 0) {
    stop(
      "sensitivityInfo() is missing required column(s): ",
      paste(missing_cols, collapse = ", ")
    )
  }

  which(info$sampleid == sample & info$treatmentid == drug)
}

pgx_extract_dose_response <- function(
  pset,
  drug,
  sample,
  summarize_replicates = TRUE
) {
  exp_idx <- pgx_resolve_experiments(pset, drug, sample)
  if (length(exp_idx) == 0) {
    stop(
      "No dose-response experiment found for sample '",
      sample,
      "' and drug '",
      drug,
      "'."
    )
  }

  raw <- pgx_sensitivity_raw(pset)
  if (is.null(dim(raw)) || length(dim(raw)) != 3 || any(dim(raw) == 0)) {
    stop(
      "No sensitivityRaw() dose-response array is available for this dataset."
    )
  }

  info <- pgx_sensitivity_info(pset)
  records <- do.call(
    rbind,
    lapply(exp_idx, function(idx) {
      dose <- as.numeric(raw[idx, , "Dose"])
      viability <- as.numeric(raw[idx, , "Viability"])
      keep <- stats::complete.cases(dose, viability)
      data.frame(
        experiment_id = rownames(info)[idx],
        sample = info$sampleid[idx],
        treatment = info$treatmentid[idx],
        dose = dose[keep],
        viability = viability[keep],
        stringsAsFactors = FALSE
      )
    })
  )

  points <- if (isTRUE(summarize_replicates)) {
    summarized <- stats::aggregate(
      viability ~ dose,
      data = records,
      FUN = median,
      na.rm = TRUE
    )
    summarized <- summarized[order(summarized$dose), , drop = FALSE]
    summarized$sample <- sample
    summarized$treatment <- drug
    summarized[, c("sample", "treatment", "dose", "viability")]
  } else {
    records[order(records$dose), , drop = FALSE]
  }

  list(
    experiment_indices = exp_idx,
    sensitivity_metadata = pgx_table_with_id(
      info[exp_idx, , drop = FALSE],
      "experiment"
    ),
    raw_points = records,
    points = points
  )
}

pgx_default_molecular_summary_stat <- function(
  pset,
  mDataType,
  summary_stat = NULL
) {
  if (!pgx_is_blank(summary_stat)) {
    return(summary_stat)
  }

  profiles <- pgx_molecular_profiles(pset)
  annotation <- tryCatch(
    S4Vectors::metadata(profiles[[mDataType]])$annotation,
    error = function(e) NA_character_
  )

  if (
    grepl("mutation|snp", annotation, ignore.case = TRUE) ||
      grepl("mutation|snp", mDataType, ignore.case = TRUE)
  ) {
    "or"
  } else {
    "mean"
  }
}

pgx_resolve_features <- function(
  pset,
  mDataType,
  features,
  match_mode = "exact",
  limit = 50
) {
  feature_info <- pgx_feature_info(pset, mDataType)
  match_fields <- colnames(feature_info)
  if (identical(pgx_resolve_match_mode(match_mode), "exact")) {
    resolution <- pgx_resolve_annotation_queries(
      feature_info,
      features,
      match_fields
    )
    matched <- resolution$rows
  } else {
    matched <- pgx_find_rows(feature_info, features, match_mode, match_fields)
    resolution <- list(
      matched_queries = pgx_clean_vector(features),
      unmatched_queries = character(),
      matched_fields = unique(unlist(strsplit(
        pgx_matched_fields(
          feature_info,
          matched,
          features,
          match_mode,
          match_fields
        ),
        ";",
        fixed = TRUE
      )))
    )
  }
  if (length(matched) == 0) {
    return(list(
      feature_ids = character(),
      matched_features = data.frame(),
      matched_queries = character(),
      unmatched_queries = pgx_clean_vector(features),
      matched_fields = character()
    ))
  }

  matched <- matched[seq_len(min(length(matched), limit))]
  feature_ids <- rownames(feature_info)[matched]
  matched_features <- pgx_table_with_id(
    feature_info[matched, , drop = FALSE],
    "feature"
  )

  list(
    feature_ids = feature_ids,
    matched_features = matched_features,
    matched_queries = resolution$matched_queries,
    unmatched_queries = resolution$unmatched_queries,
    matched_fields = resolution$matched_fields
  )
}

pgx_get_molecular_matrix <- function(
  pset,
  mDataType,
  features,
  samples = NULL,
  feature_match_mode = "exact",
  summary_stat = NULL,
  limit_features = 25,
  limit_samples = 25
) {
  available_profiles <- pgx_mdata_names(pset)
  if (!mDataType %in% available_profiles) {
    available_text <- if (length(available_profiles) == 0) {
      "none"
    } else {
      paste(available_profiles, collapse = ", ")
    }
    stop(
      "Molecular profile '",
      mDataType,
      "' is not available. Available profiles are: ",
      available_text,
      ". Call pgx_list_available_covariates() before requesting molecular data."
    )
  }

  samples <- pgx_clean_vector(samples)
  if (is.null(samples)) {
    samples <- pgx_sample_names(pset)
  }
  samples <- samples[seq_len(min(
    length(samples),
    max(1L, as.integer(limit_samples))
  ))]
  summary_stat <- pgx_default_molecular_summary_stat(
    pset,
    mDataType,
    summary_stat
  )

  resolved <- pgx_resolve_features(
    pset = pset,
    mDataType = mDataType,
    features = features,
    match_mode = feature_match_mode,
    limit = max(1L, as.integer(limit_features))
  )
  if (length(resolved$feature_ids) == 0) {
    return(list(
      matrix = matrix(numeric(), nrow = 0, ncol = 0),
      matched_features = data.frame(),
      summary_stat = summary_stat
    ))
  }

  mat <- pgx_summarize_molecular_profiles(
    pset,
    mDataType = mDataType,
    features = resolved$feature_ids,
    cell.lines = samples,
    summary.stat = summary_stat
  )

  list(
    matrix = mat,
    matched_features = resolved$matched_features,
    summary_stat = summary_stat
  )
}

pgx_matrix_records <- function(mat, row_name = "feature", col_name = "sample") {
  if (length(dim(mat)) != 2) {
    stop("Expected a two-dimensional molecular profile matrix.")
  }
  if (any(dim(mat) == 0)) {
    records <- data.frame(
      row = character(),
      column = character(),
      value = numeric(),
      stringsAsFactors = FALSE
    )
    names(records) <- c(row_name, col_name, "value")
    return(records)
  }

  records <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
  names(records) <- c(row_name, col_name, "value")
  records
}

pgx_required_column_result <- function(table_name, provided, required) {
  provided <- pgx_clean_vector(provided)
  if (is.null(provided)) {
    provided <- character()
  }
  missing <- setdiff(required, provided)
  data.frame(
    table = table_name,
    required = paste(required, collapse = ", "),
    provided_count = length(provided),
    missing = paste(missing, collapse = ", "),
    ok = length(missing) == 0,
    stringsAsFactors = FALSE
  )
}

pgx_pair_table <- function(pset) {
  info <- pgx_sensitivity_info(pset)
  if (all(c("sampleid", "treatmentid") %in% colnames(info))) {
    return(unique(data.frame(
      sample = info$sampleid,
      treatment = info$treatmentid,
      stringsAsFactors = FALSE
    )))
  }

  if (all(c("sampleid", "treatment1id", "treatment2id") %in% colnames(info))) {
    return(unique(data.frame(
      sample = info$sampleid,
      treatment = paste(info$treatment1id, info$treatment2id, sep = " + "),
      stringsAsFactors = FALSE
    )))
  }

  treatment_response <- tryCatch(
    pgx_treatment_response(pset),
    error = function(e) {
      list()
    }
  )
  for (table_name in names(treatment_response)) {
    x <- as.data.frame(treatment_response[[table_name]])
    if (all(c("sampleid", "treatment1id", "treatment2id") %in% colnames(x))) {
      return(unique(data.frame(
        sample = x$sampleid,
        treatment = paste(x$treatment1id, x$treatment2id, sep = " + "),
        stringsAsFactors = FALSE
      )))
    }
  }

  data.frame(sample = character(), treatment = character())
}

pgx_sensitivity_measures <- function(pset) {
  response <- tryCatch(
    pgx_treatment_response(pset),
    error = function(e) NULL
  )
  if (!is.null(response) && "mono_profiles" %in% names(response)) {
    return(pgx_numeric_metric_columns(response$mono_profiles))
  }
  if (pgx_has_legacy_attrs(pset) && !is.null(response$profiles)) {
    return(pgx_numeric_metric_columns(response$profiles))
  }

  measures <- tryCatch(
    suppressWarnings(PharmacoGx::sensitivityMeasures(pset)),
    error = function(e) NULL
  )
  if (!is.null(measures)) {
    return(measures)
  }

  if (is.null(response)) {
    response <- pgx_treatment_response(pset)
  }
  profiles <- response$profiles
  if (is.null(profiles)) {
    return(character())
  }

  pgx_numeric_metric_columns(profiles)
}

pgx_sensitivity_info <- function(pset) {
  if (pgx_has_legacy_attrs(pset)) {
    sensitivity <- pgx_legacy_attr(pset, c("treatmentResponse", "sensitivity"))
    if (!is.null(sensitivity) && "info" %in% names(sensitivity)) {
      return(pgx_as_data_frame(sensitivity[["info"]]))
    }
  }

  info <- tryCatch(
    suppressWarnings(PharmacoGx::sensitivityInfo(pset)),
    error = function(e) NULL
  )
  if (!is.null(info)) {
    return(pgx_as_data_frame(info))
  }

  sensitivity <- pgx_legacy_attr(pset, c("treatmentResponse", "sensitivity"))
  if (is.null(sensitivity)) {
    return(data.frame())
  }
  candidates <- intersect(
    c("info", "assay_metadata", "response", "profiles"),
    names(sensitivity)
  )
  if (length(candidates) == 0) {
    return(data.frame())
  }
  pgx_as_data_frame(sensitivity[[candidates[[1]]]])
}

pgx_sensitivity_raw <- function(pset) {
  if (pgx_has_legacy_attrs(pset)) {
    sensitivity <- pgx_legacy_attr(pset, c("treatmentResponse", "sensitivity"))
    return(sensitivity$raw)
  }

  raw <- tryCatch(
    suppressWarnings(PharmacoGx::sensitivityRaw(pset)),
    error = function(e) NULL
  )
  if (!is.null(raw)) {
    return(raw)
  }

  sensitivity <- pgx_legacy_attr(pset, c("treatmentResponse", "sensitivity"))
  sensitivity$raw
}

pgx_load_many_datasets <- function(datasets) {
  datasets <- pgx_clean_vector(datasets)
  if (is.null(datasets) || length(datasets) < 2) {
    stop("At least two datasets are required.")
  }
  stats::setNames(lapply(datasets, pgx_load_dataset), datasets)
}

pgx_treatment_response <- function(pset) {
  if (pgx_has_legacy_attrs(pset)) {
    response <- pgx_legacy_attr(pset, c("treatmentResponse", "sensitivity"))
    if (is.null(response)) {
      stop("Could not access legacy sensitivity data.")
    }
    return(response[vapply(
      response,
      function(x) {
        is.data.frame(x) || is.matrix(x) || length(dim(x)) >= 2
      },
      logical(1)
    )])
  }

  response <- tryCatch(
    suppressWarnings(PharmacoGx::treatmentResponse(pset)),
    error = function(e) {
      NULL
    }
  )
  if (!is.null(response)) {
    return(response)
  }

  response <- pgx_legacy_attr(pset, c("treatmentResponse", "sensitivity"))
  if (is.null(response)) {
    stop("Could not access treatmentResponse() or legacy sensitivity data.")
  }

  response[vapply(
    response,
    function(x) {
      is.data.frame(x) || is.matrix(x) || length(dim(x)) >= 2
    },
    logical(1)
  )]
}

pgx_treatment_response_table_names <- function(pset) {
  names(pgx_treatment_response(pset))
}

pgx_response_table <- function(
  pset,
  table = "auto",
  preferred = c("profiles", "synergy", "raw"),
  as_data_table = FALSE
) {
  treatment_response <- pgx_treatment_response(pset)
  tables <- names(treatment_response)
  if (length(tables) == 0) {
    stop("No treatmentResponse() tables are available.")
  }

  table <- if (pgx_is_blank(table) || identical(table, "auto")) {
    preferred_match <- intersect(preferred, tables)
    if (length(preferred_match) > 0) {
      preferred_match[[1]]
    } else {
      tables[[1]]
    }
  } else {
    table
  }

  if (!table %in% tables) {
    stop(
      "Unknown treatmentResponse table '",
      table,
      "'. Available tables are: ",
      paste(tables, collapse = ", ")
    )
  }

  x <- treatment_response[[table]]
  if (isTRUE(as_data_table)) {
    if (!requireNamespace("data.table", quietly = TRUE)) {
      stop("The data.table package is required for combo table summaries.")
    }
    x <- data.table::as.data.table(x)
  } else {
    x <- as.data.frame(x)
  }

  list(name = table, table = x)
}

pgx_score_candidates <- function(x) {
  cols <- colnames(x)
  numeric_cols <- cols[vapply(x, is.numeric, logical(1))]
  candidates <- grep(
    paste(
      c(
        "synergy",
        "score",
        "delta",
        "excess",
        "HSA",
        "Bliss",
        "Loewe",
        "ZIP",
        "combo_aac",
        "aac",
        "viability"
      ),
      collapse = "|"
    ),
    numeric_cols,
    value = TRUE,
    ignore.case = TRUE
  )
  candidates <- candidates[
    !grepl("missing|expected", candidates, ignore.case = TRUE)
  ]

  priority <- c(
    "ZIP_delta",
    "ZIP_score",
    "Bliss_score",
    "HSA_score",
    "Loewe_score",
    "mean_HSA_score_recomputed",
    "median_HSA_score_recomputed",
    "mean_HSA_score_obs",
    "median_HSA_score_obs",
    "mean_HSA_score",
    "median_HSA_score",
    "HSA_synergy_score",
    "HSA_excess",
    "mean_Bliss_score_recomputed",
    "median_Bliss_score_recomputed",
    "mean_Bliss_score_obs",
    "median_Bliss_score_obs",
    "mean_Bliss_score",
    "median_Bliss_score",
    "Bliss_synergy_score",
    "Bliss_excess",
    "mean_Loewe_score",
    "median_Loewe_score",
    "mean_ZIP_score",
    "median_ZIP_score",
    "ZIP_synergy_score",
    "combo_aac_observed"
  )

  unique(c(intersect(priority, candidates), candidates))
}

pgx_choose_score_column <- function(x, score = NULL) {
  candidates <- pgx_score_candidates(x)
  if (!pgx_is_blank(score)) {
    if (!score %in% colnames(x)) {
      stop(
        "Unknown score column '",
        score,
        "'. Available columns are: ",
        paste(colnames(x), collapse = ", ")
      )
    }
    if (!is.numeric(x[[score]])) {
      stop("Score column '", score, "' must be numeric.")
    }
    return(score)
  }

  if (length(candidates) == 0) {
    stop("No numeric synergy or response score columns were detected.")
  }

  candidates[[1]]
}

pgx_score_direction <- function(score, rank_direction = "auto") {
  rank_direction <- match.arg(rank_direction, c("auto", "lowest", "highest"))
  if (!identical(rank_direction, "auto")) {
    return(rank_direction)
  }

  if (grepl("viability|ic50|ec50|rmse", score, ignore.case = TRUE)) {
    "lowest"
  } else {
    "highest"
  }
}

pgx_filter_combo_rows <- function(
  x,
  treatment1 = NULL,
  treatment2 = NULL,
  match_mode = "exact",
  symmetric = TRUE
) {
  required <- c("treatment1id", "treatment2id")
  missing <- setdiff(required, colnames(x))
  if (length(missing) > 0) {
    stop(
      "Combo table is missing required column(s): ",
      paste(missing, collapse = ", ")
    )
  }

  treatment1 <- pgx_clean_vector(treatment1)
  treatment2 <- pgx_clean_vector(treatment2)
  if (is.null(treatment1) && is.null(treatment2)) {
    return(x)
  }

  match_mode <- pgx_resolve_match_mode(match_mode)
  t1 <- as.character(x$treatment1id)
  t2 <- as.character(x$treatment2id)

  if (!is.null(treatment1) && !is.null(treatment2) && isTRUE(symmetric)) {
    keep <- (pgx_match_filter_values(t1, treatment1, match_mode) &
      pgx_match_filter_values(t2, treatment2, match_mode)) |
      (pgx_match_filter_values(t1, treatment2, match_mode) &
        pgx_match_filter_values(t2, treatment1, match_mode))
  } else {
    keep <- rep(TRUE, nrow(x))
    if (!is.null(treatment1)) {
      keep <- keep & pgx_match_filter_values(t1, treatment1, match_mode)
    }
    if (!is.null(treatment2)) {
      keep <- keep & pgx_match_filter_values(t2, treatment2, match_mode)
    }
  }

  x[keep, , drop = FALSE]
}

pgx_feature_response_association <- function(
  mat,
  response_values,
  method,
  p_adjust_method
) {
  if (nrow(mat) == 0 || ncol(mat) == 0) {
    return(data.frame())
  }

  common_samples <- intersect(names(response_values), colnames(mat))

  results <- do.call(
    rbind,
    lapply(rownames(mat), function(feature_id) {
      x <- if (length(common_samples) > 0) {
        as.vector(mat[feature_id, common_samples, drop = TRUE])
      } else {
        vector()
      }
      y <- as.numeric(response_values[common_samples])
      complete <- stats::complete.cases(x, y) & is.finite(y)
      if (is.numeric(x)) {
        complete <- complete & is.finite(x)
      }
      x <- x[complete]
      y <- y[complete]

      if (length(y) < 3 || length(unique(x)) < 2 || length(unique(y)) < 2) {
        return(data.frame(
          feature = feature_id,
          method = method,
          n = length(y),
          effect_size = NA_real_,
          p_value = NA_real_,
          stringsAsFactors = FALSE
        ))
      }

      if (method %in% c("spearman", "pearson") && is.numeric(x)) {
        test <- stats::cor.test(x, y, method = method, exact = FALSE)
        effect <- unname(test$estimate)
        p_value <- test$p.value
        method_used <- method
      } else {
        group <- as.factor(x)
        if (length(levels(group)) != 2) {
          return(data.frame(
            feature = feature_id,
            method = "wilcoxon",
            n = length(y),
            effect_size = NA_real_,
            p_value = NA_real_,
            stringsAsFactors = FALSE
          ))
        }
        test <- suppressWarnings(stats::wilcox.test(y ~ group))
        effect <- unname(diff(tapply(y, group, median, na.rm = TRUE))[[1]])
        p_value <- test$p.value
        method_used <- "wilcoxon"
      }

      data.frame(
        feature = feature_id,
        method = method_used,
        n = length(y),
        effect_size = effect,
        p_value = p_value,
        stringsAsFactors = FALSE
      )
    })
  )

  results$fdr <- stats::p.adjust(results$p_value, method = p_adjust_method)
  results$p_value_method <- ifelse(
    results$method == "spearman",
    "asymptotic (exact=FALSE)",
    ifelse(results$method == "pearson", "Pearson t-test", "wilcoxon default")
  )
  results
}

pgx_fisher_p <- function(p_values) {
  p_values <- p_values[!is.na(p_values) & is.finite(p_values)]
  if (length(p_values) == 0) {
    return(NA_real_)
  }

  p_values <- pmax(p_values, .Machine$double.xmin)
  stats::pchisq(
    -2 * sum(log(p_values)),
    df = 2 * length(p_values),
    lower.tail = FALSE
  )
}

pgx_known_biomarker_candidates <- data.frame(
  drug = c("Gemcitabine", "Gemcitabine", "Anagrelide"),
  drug_aliases = c("gemcitabine;dFdC", "gemcitabine;dFdC", "anagrelide"),
  feature = c("SLC29A1", "SLC28A1", "PDE3A"),
  feature_aliases = c("hENT1;ENT1", "hCNT1;CNT1", "PDE3A"),
  mDataType_hint = c("rna", "rna", "rna;mutation;cnv"),
  context = c(
    "Nucleoside transporter candidate for gemcitabine uptake and response.",
    "Nucleoside transporter candidate sometimes evaluated with gemcitabine response.",
    "Mechanistic target candidate for anagrelide response."
  ),
  stringsAsFactors = FALSE
)

pgx_capture <- function(expr) {
  warnings <- character()
  messages <- character()

  value <- withCallingHandlers(
    tryCatch(
      expr,
      error = function(e) {
        structure(
          list(message = conditionMessage(e)),
          class = "pgx_tool_error"
        )
      }
    ),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      messages <<- c(messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )

  list(
    ok = !inherits(value, "pgx_tool_error"),
    value = value,
    warnings = unique(warnings),
    messages = unique(messages)
  )
}

pgx_metric_value <- function(expr) {
  result <- pgx_capture(expr)

  if (!isTRUE(result$ok)) {
    return(list(value = NA_real_, error = result$value$message))
  }

  list(
    value = unname(as.numeric(result$value)[1]),
    warnings = result$warnings,
    messages = result$messages
  )
}

pgx_requirement_status <- function(requested, matched) {
  if (pgx_is_blank(requested)) {
    return(NA)
  }
  length(matched) > 0
}

pgx_unloadable_requirements_result <- function(
  dataset,
  drug = NULL,
  sample_query = NULL,
  metric = NULL,
  mDataType = NULL,
  features = NULL,
  error,
  check_mode
) {
  requested_names <- c(
    if (!pgx_is_blank(drug)) "drug",
    if (!pgx_is_blank(sample_query)) "sample_context",
    if (!pgx_is_blank(metric)) "response_metric",
    if (!pgx_is_blank(mDataType)) "molecular_profile",
    if (!pgx_is_blank(features)) "molecular_feature"
  )
  data.frame(
    dataset = dataset,
    check_mode = check_mode,
    loadable = FALSE,
    drug_available = NA,
    drug_matches = "",
    drug_match_fields = "",
    drug_unmatched_queries = "",
    sample_available = NA,
    sample_match_count = 0L,
    sample_match_fields = "",
    sample_query_variants = "",
    metric_available = NA,
    metric_matches = "",
    metric_match_type = "",
    molecular_profile_available = NA,
    molecular_profile_usable = NA,
    molecular_profile_matches = "",
    feature_available = NA,
    feature_matches = "",
    feature_matched_queries = "",
    feature_unmatched_queries = "",
    requirements_requested = length(requested_names),
    requirements_satisfied = 0L,
    requirements_fraction = 0,
    missing_requirements = paste(
      c("not_loadable", requested_names),
      collapse = ";"
    ),
    requirements_met = FALSE,
    error = error,
    stringsAsFactors = FALSE
  )
}

pgx_check_loaded_dataset_requirements <- function(
  dataset,
  drug = NULL,
  sample_query = NULL,
  metric = NULL,
  mDataType = NULL,
  features = NULL
) {
  loaded <- pgx_capture(pgx_load_dataset(dataset))
  if (!isTRUE(loaded$ok)) {
    return(pgx_unloadable_requirements_result(
      dataset = dataset,
      drug = drug,
      sample_query = sample_query,
      metric = metric,
      mDataType = mDataType,
      features = features,
      error = loaded$value$message,
      check_mode = "in_process"
    ))
  }

  pset <- loaded$value
  drug_resolution <- if (pgx_is_blank(drug)) {
    list(
      identifiers = character(),
      matched_fields = character(),
      unmatched_queries = character()
    )
  } else {
    tryCatch(
      pgx_resolve_treatment_identifiers(pset, drug),
      error = function(e) {
        list(
          identifiers = character(),
          matched_fields = character(),
          unmatched_queries = pgx_clean_vector(drug)
        )
      }
    )
  }
  drug_matches <- drug_resolution$identifiers

  sample_info <- tryCatch(pgx_sample_info(pset), error = function(e) {
    data.frame()
  })
  tissue_fields <- grep(
    "tissue|site|histology|lineage|disease|cancer|origin|panel",
    colnames(sample_info),
    ignore.case = TRUE,
    value = TRUE
  )
  sample_variants <- pgx_sample_context_variants(sample_query)
  sample_rows <- if (pgx_is_blank(sample_query) || length(tissue_fields) == 0) {
    integer()
  } else {
    pgx_find_rows(sample_info, sample_variants, "contains", tissue_fields)
  }
  sample_match_fields <- if (length(sample_rows) == 0) {
    character()
  } else {
    unique(unlist(strsplit(
      pgx_matched_fields(
        sample_info,
        sample_rows,
        sample_variants,
        "contains",
        tissue_fields
      ),
      ";",
      fixed = TRUE
    )))
  }

  available_metrics <- tryCatch(
    pgx_sensitivity_measures(pset),
    error = function(e) character()
  )
  metric_resolution <- pgx_resolve_metric(metric, available_metrics)
  metric_match <- if (pgx_is_blank(metric)) {
    character()
  } else if (!is.null(metric_resolution$resolved)) {
    metric_resolution$resolved
  } else {
    character()
  }

  capabilities <- tryCatch(
    pgx_molecular_profile_capabilities(pset),
    error = function(e) data.frame()
  )
  profile_row <- pgx_resolve_molecular_profile(mDataType, capabilities)
  profile_available <- nrow(profile_row) > 0
  usable_profile_rows <- profile_row[profile_row$usable, , drop = FALSE]
  profile_usable <- nrow(usable_profile_rows) > 0
  selected_profile <- if (profile_usable) {
    usable_profile_rows$molecular_profile[[1]]
  } else if (profile_available) {
    profile_row$molecular_profile[[1]]
  } else {
    NULL
  }
  feature_matches <- character()
  feature_matched_queries <- character()
  feature_unmatched_queries <- pgx_clean_vector(features)
  if (!pgx_is_blank(features) && profile_usable) {
    feature_resolution <- tryCatch(
      pgx_resolve_features(
        pset,
        mDataType = selected_profile,
        features = features,
        match_mode = "exact",
        limit = 100
      ),
      error = function(e) {
        list(
          feature_ids = character(),
          matched_queries = character(),
          unmatched_queries = pgx_clean_vector(features)
        )
      }
    )
    feature_matches <- feature_resolution$feature_ids
    feature_matched_queries <- feature_resolution$matched_queries
    feature_unmatched_queries <- feature_resolution$unmatched_queries
  }

  requested_statuses <- c(
    drug = pgx_requirement_status(drug, drug_matches),
    sample_context = pgx_requirement_status(sample_query, sample_rows),
    response_metric = pgx_requirement_status(metric, metric_match),
    molecular_profile = if (pgx_is_blank(mDataType)) NA else profile_usable,
    molecular_feature = pgx_requirement_status(features, feature_matches)
  )
  requested <- requested_statuses[!is.na(requested_statuses)]
  requirements_requested <- length(requested)
  requirements_satisfied <- sum(requested)
  requirements_fraction <- if (requirements_requested == 0) {
    1
  } else {
    requirements_satisfied / requirements_requested
  }
  data.frame(
    dataset = dataset,
    check_mode = "in_process",
    loadable = TRUE,
    drug_available = pgx_requirement_status(drug, drug_matches),
    drug_matches = paste(head(drug_matches, 10), collapse = ";"),
    drug_match_fields = paste(drug_resolution$matched_fields, collapse = ";"),
    drug_unmatched_queries = paste(
      drug_resolution$unmatched_queries,
      collapse = ";"
    ),
    sample_available = pgx_requirement_status(sample_query, sample_rows),
    sample_match_count = length(sample_rows),
    sample_match_fields = paste(sample_match_fields, collapse = ";"),
    sample_query_variants = paste(sample_variants, collapse = ";"),
    metric_available = pgx_requirement_status(metric, metric_match),
    metric_matches = paste(metric_match, collapse = ";"),
    metric_match_type = if (pgx_is_blank(metric)) {
      ""
    } else if (!is.null(metric_resolution$match_type)) {
      metric_resolution$match_type
    } else {
      "unresolved"
    },
    molecular_profile_available = if (pgx_is_blank(mDataType)) {
      NA
    } else {
      profile_available
    },
    molecular_profile_usable = if (pgx_is_blank(mDataType)) {
      NA
    } else {
      profile_usable
    },
    molecular_profile_matches = paste(
      profile_row$molecular_profile,
      collapse = ";"
    ),
    feature_available = pgx_requirement_status(features, feature_matches),
    feature_matches = paste(head(feature_matches, 10), collapse = ";"),
    feature_matched_queries = paste(feature_matched_queries, collapse = ";"),
    feature_unmatched_queries = paste(
      feature_unmatched_queries,
      collapse = ";"
    ),
    requirements_requested = requirements_requested,
    requirements_satisfied = requirements_satisfied,
    requirements_fraction = requirements_fraction,
    missing_requirements = paste(names(requested)[!requested], collapse = ";"),
    requirements_met = all(requested),
    error = "",
    stringsAsFactors = FALSE
  )
}

pgx_check_capability_snapshot_requirements <- function(
  snapshot,
  dataset,
  drug = NULL,
  sample_query = NULL,
  metric = NULL,
  mDataType = NULL,
  features = NULL,
  cache_status = "generated"
) {
  treatment_info <- snapshot$treatment_info
  drug_resolution <- if (pgx_is_blank(drug)) {
    list(
      rows = integer(),
      matched_fields = character(),
      matched_queries = character(),
      unmatched_queries = character()
    )
  } else {
    pgx_resolve_annotation_queries(
      treatment_info,
      drug,
      search_fields = colnames(treatment_info)
    )
  }
  drug_matches <- pgx_first_row_identifier(
    treatment_info,
    drug_resolution$rows,
    c("treatmentid", "drugid", "treatment", "drug", "name", "compound")
  )

  sample_info <- snapshot$sample_info
  tissue_fields <- grep(
    "tissue|site|histology|lineage|disease|cancer|origin|panel",
    colnames(sample_info),
    ignore.case = TRUE,
    value = TRUE
  )
  sample_variants <- pgx_sample_context_variants(sample_query)
  sample_rows <- if (pgx_is_blank(sample_query) || length(tissue_fields) == 0) {
    integer()
  } else {
    pgx_find_rows(sample_info, sample_variants, "contains", tissue_fields)
  }
  sample_match_fields <- if (length(sample_rows) == 0) {
    character()
  } else {
    unique(unlist(strsplit(
      pgx_matched_fields(
        sample_info,
        sample_rows,
        sample_variants,
        "contains",
        tissue_fields
      ),
      ";",
      fixed = TRUE
    )))
  }

  available_metrics <- snapshot$sensitivity_measures
  metric_resolution <- pgx_resolve_metric(metric, available_metrics)
  metric_match <- if (pgx_is_blank(metric)) {
    character()
  } else if (!is.null(metric_resolution$resolved)) {
    metric_resolution$resolved
  } else {
    character()
  }

  capabilities <- snapshot$molecular_capabilities
  profile_row <- pgx_resolve_molecular_profile(mDataType, capabilities)
  profile_available <- nrow(profile_row) > 0
  usable_profile_rows <- profile_row[profile_row$usable, , drop = FALSE]
  profile_usable <- nrow(usable_profile_rows) > 0
  selected_profile <- if (profile_usable) {
    usable_profile_rows$molecular_profile[[1]]
  } else if (profile_available) {
    profile_row$molecular_profile[[1]]
  } else {
    NULL
  }

  feature_matches <- character()
  feature_matched_queries <- character()
  feature_unmatched_queries <- pgx_clean_vector(features)
  if (!pgx_is_blank(features) && profile_usable) {
    feature_info <- snapshot$feature_info[[selected_profile]]
    feature_resolution <- pgx_resolve_annotation_queries(
      feature_info,
      features,
      search_fields = colnames(feature_info)
    )
    feature_matches <- pgx_first_row_identifier(
      feature_info,
      feature_resolution$rows,
      c(
        "featureid",
        "geneid",
        "gene_id",
        "symbol",
        "gene_symbol",
        "ensembl",
        "probe",
        "name"
      )
    )
    feature_matched_queries <- feature_resolution$matched_queries
    feature_unmatched_queries <- feature_resolution$unmatched_queries
  }

  requested_statuses <- c(
    drug = pgx_requirement_status(drug, drug_matches),
    sample_context = pgx_requirement_status(sample_query, sample_rows),
    response_metric = pgx_requirement_status(metric, metric_match),
    molecular_profile = if (pgx_is_blank(mDataType)) NA else profile_usable,
    molecular_feature = pgx_requirement_status(features, feature_matches)
  )
  requested <- requested_statuses[!is.na(requested_statuses)]
  requirements_requested <- length(requested)
  requirements_satisfied <- sum(requested)
  requirements_fraction <- if (requirements_requested == 0) {
    1
  } else {
    requirements_satisfied / requirements_requested
  }

  data.frame(
    dataset = dataset,
    check_mode = paste0("isolated_worker_", cache_status),
    loadable = TRUE,
    drug_available = pgx_requirement_status(drug, drug_matches),
    drug_matches = paste(head(drug_matches, 10), collapse = ";"),
    drug_match_fields = paste(drug_resolution$matched_fields, collapse = ";"),
    drug_unmatched_queries = paste(
      drug_resolution$unmatched_queries,
      collapse = ";"
    ),
    sample_available = pgx_requirement_status(sample_query, sample_rows),
    sample_match_count = length(sample_rows),
    sample_match_fields = paste(sample_match_fields, collapse = ";"),
    sample_query_variants = paste(sample_variants, collapse = ";"),
    metric_available = pgx_requirement_status(metric, metric_match),
    metric_matches = paste(metric_match, collapse = ";"),
    metric_match_type = if (pgx_is_blank(metric)) {
      ""
    } else if (!is.null(metric_resolution$match_type)) {
      metric_resolution$match_type
    } else {
      "unresolved"
    },
    molecular_profile_available = if (pgx_is_blank(mDataType)) {
      NA
    } else {
      profile_available
    },
    molecular_profile_usable = if (pgx_is_blank(mDataType)) {
      NA
    } else {
      profile_usable
    },
    molecular_profile_matches = paste(
      profile_row$molecular_profile,
      collapse = ";"
    ),
    feature_available = pgx_requirement_status(features, feature_matches),
    feature_matches = paste(head(feature_matches, 10), collapse = ";"),
    feature_matched_queries = paste(feature_matched_queries, collapse = ";"),
    feature_unmatched_queries = paste(
      feature_unmatched_queries,
      collapse = ";"
    ),
    requirements_requested = requirements_requested,
    requirements_satisfied = requirements_satisfied,
    requirements_fraction = requirements_fraction,
    missing_requirements = paste(names(requested)[!requested], collapse = ";"),
    requirements_met = all(requested),
    error = "",
    stringsAsFactors = FALSE
  )
}

pgx_check_dataset_requirements <- function(
  dataset,
  drug = NULL,
  sample_query = NULL,
  metric = NULL,
  mDataType = NULL,
  features = NULL
) {
  local_path <- pgx_local_dataset_path(dataset)
  if (is.null(local_path) || pgx_env_flag("PHARMACOGX_MCP_IN_PROBE_WORKER")) {
    return(pgx_check_loaded_dataset_requirements(
      dataset = dataset,
      drug = drug,
      sample_query = sample_query,
      metric = metric,
      mDataType = mDataType,
      features = features
    ))
  }

  probed <- pgx_probe_local_pset_isolated(local_path)
  if (!isTRUE(probed$ok)) {
    return(pgx_unloadable_requirements_result(
      dataset = dataset,
      drug = drug,
      sample_query = sample_query,
      metric = metric,
      mDataType = mDataType,
      features = features,
      error = probed$message,
      check_mode = "isolated_worker_failed"
    ))
  }
  pgx_check_capability_snapshot_requirements(
    snapshot = probed$snapshot,
    dataset = dataset,
    drug = drug,
    sample_query = sample_query,
    metric = metric,
    mDataType = mDataType,
    features = features,
    cache_status = probed$cache_status
  )
}

pgx_check_analysis_requirements <- function(
  datasets = NULL,
  drug = NULL,
  sample_query = NULL,
  metric = NULL,
  mDataType = NULL,
  features = NULL,
  include_local = TRUE,
  include_bundled = TRUE,
  limit = 20
) {
  datasets <- pgx_clean_vector(datasets)
  if (is.null(datasets)) {
    local <- if (isTRUE(include_local)) {
      pgx_local_pset_catalog(include_missing = FALSE)$name
    } else {
      character()
    }
    bundled <- if (isTRUE(include_bundled)) {
      pgx_supported_datasets
    } else {
      character()
    }
    datasets <- unique(c(local, bundled))
  }
  limit <- max(1L, as.integer(limit))
  datasets <- head(datasets, limit)
  checks <- do.call(
    rbind,
    lapply(
      datasets,
      pgx_check_dataset_requirements,
      drug = drug,
      sample_query = sample_query,
      metric = metric,
      mDataType = mDataType,
      features = features
    )
  )
  matching <- checks$dataset[checks$requirements_met]
  ranked <- checks[
    order(
      checks$loadable,
      checks$requirements_fraction,
      checks$requirements_satisfied,
      decreasing = TRUE
    ),
    ,
    drop = FALSE
  ]
  candidate_fields <- c(
    "dataset",
    "check_mode",
    "loadable",
    "requirements_satisfied",
    "requirements_requested",
    "requirements_fraction",
    "missing_requirements",
    "drug_matches",
    "sample_match_count",
    "metric_matches",
    "molecular_profile_matches",
    "feature_matches",
    "error"
  )
  best_candidates <- head(
    ranked[ranked$loadable, candidate_fields, drop = FALSE],
    5L
  )
  unchecked <- checks[!checks$loadable, c("dataset", "error"), drop = FALSE]

  list(
    requested = list(
      drug = drug,
      sample_query = sample_query,
      metric = metric,
      molecular_profile = mDataType,
      features = pgx_clean_vector(features)
    ),
    checked = nrow(checks),
    compatible_datasets = matching,
    best_local_candidates = best_candidates,
    unchecked_datasets = unchecked,
    requirements = checks,
    local_capability_checks_are_isolated = TRUE,
    remote_catalog_assesses_compatibility = FALSE,
    next_step = if (length(matching) > 0) {
      paste(
        "Use one compatible dataset with pgx_list_available_covariates(),",
        "then run the requested analysis."
      )
    } else {
      paste(
        "No successfully checked dataset satisfies every requirement. Review",
        "best_local_candidates and the per-requirement identifier matches,",
        "then inspect covariates for the closest candidates before ruling them",
        "out. Entries in unchecked_datasets were not evaluated and may still",
        "be relevant. The downloadable catalog does not report enough",
        "treatment, sample, endpoint, and molecular coverage to establish",
        "compatibility; do not use it as the automatic next step."
      )
    }
  )
}

pgx_list_example_datasets <- function() {
  list(
    datasets = pgx_dataset_catalog(),
    note = paste(
      "Bundled datasets are toy fixtures. Local manifest aliases point to",
      "real PSet-backed demo objects when those files exist. Use downloadPSet()",
      "only when an explicit large external dataset workflow is intended."
    )
  )
}

pgx_list_local_psets <- function(
  manifest_path = NULL,
  roots = NULL,
  discover = TRUE,
  include_missing = TRUE
) {
  manifest_path <- if (pgx_is_blank(manifest_path)) {
    pgx_default_local_manifest_path()
  } else {
    path.expand(manifest_path)
  }
  catalog <- pgx_local_pset_catalog(
    manifest_path = manifest_path,
    roots = roots,
    discover = discover,
    include_missing = include_missing
  )
  discovery_roots <- if (pgx_is_blank(roots)) {
    pgx_default_local_pset_roots()
  } else {
    normalizePath(path.expand(pgx_clean_vector(roots)), mustWork = FALSE)
  }

  list(
    manifest_path = normalizePath(manifest_path, mustWork = FALSE),
    discovery_roots = discovery_roots,
    automatic_discovery = isTRUE(discover),
    include_missing = isTRUE(include_missing),
    total = nrow(catalog),
    local_psets = catalog,
    note = paste(
      "Use the name column as the dataset argument. Configure manifest entries",
      "with PHARMACOGX_LOCAL_PSET_MANIFEST or discovery roots with",
      "PHARMACOGX_LOCAL_PSET_DIRS."
    )
  )
}

pgx_list_entities <- function(
  dataset = "GDSCsmall",
  entity = "treatments",
  query = NULL,
  match_mode = "contains",
  offset = 0,
  limit = 25
) {
  entity <- match.arg(
    entity,
    c("treatments", "samples", "sensitivity_measures")
  )
  match_mode <- pgx_resolve_match_mode(match_mode)
  offset <- max(0L, as.integer(offset))
  limit <- max(1L, as.integer(limit))

  pset <- pgx_load_dataset(dataset)
  values <- switch(
    entity,
    "treatments" = pgx_treatment_names(pset),
    "samples" = pgx_sample_names(pset),
    "sensitivity_measures" = pgx_sensitivity_measures(pset)
  )
  values <- sort(unique(pgx_normalize_text(values)))
  query <- pgx_clean_vector(query)
  if (!is.null(query)) {
    keep <- Reduce(
      `|`,
      lapply(query, function(value) {
        if (identical(match_mode, "exact")) {
          pgx_casefold(values) == pgx_casefold(value)
        } else {
          grepl(pgx_casefold(value), pgx_casefold(values), fixed = TRUE)
        }
      })
    )
    values <- values[keep]
  }
  selected <- if (offset >= length(values)) {
    character()
  } else {
    start <- offset + 1L
    values[seq.int(start, min(length(values), offset + limit))]
  }

  list(
    dataset = dataset,
    entity = entity,
    query = query,
    match_mode = match_mode,
    offset = offset,
    total = length(values),
    returned = length(selected),
    next_offset = if (offset + length(selected) < length(values)) {
      offset + length(selected)
    } else {
      NULL
    },
    values = selected
  )
}

pgx_list_available_covariates <- function(dataset = "GDSCsmall") {
  pset <- pgx_load_dataset(dataset)
  sensitivity_info <- pgx_sensitivity_info(pset)
  treatment_response_tables <- tryCatch(
    pgx_treatment_response_table_names(pset),
    error = function(e) character()
  )
  molecular_profiles <- tryCatch(
    pgx_mdata_names(pset),
    error = function(e) character()
  )

  feature_fields <- lapply(molecular_profiles, function(mDataType) {
    tryCatch(
      colnames(pgx_feature_info(pset, mDataType)),
      error = function(e) character()
    )
  })
  names(feature_fields) <- molecular_profiles
  molecular_capabilities <- pgx_molecular_profile_capabilities(pset)

  list(
    dataset = dataset,
    sample_fields = colnames(pgx_sample_info(pset)),
    treatment_fields = colnames(pgx_treatment_info(pset)),
    sensitivity_info_fields = colnames(sensitivity_info),
    sensitivity_measures = pgx_sensitivity_measures(pset),
    treatment_response_tables = treatment_response_tables,
    molecular_profiles = molecular_profiles,
    molecular_profile_capabilities = molecular_capabilities,
    feature_fields = feature_fields,
    note = paste(
      "Use these field names with metadata and filter tools. Field names differ",
      "between PharmacoSets, so discover covariates before filtering."
    )
  )
}

pgx_get_sample_metadata <- function(
  dataset = "GDSCsmall",
  samples,
  fields = NULL,
  match_mode = "exact",
  limit = 25
) {
  pset <- pgx_load_dataset(dataset)
  info <- pgx_sample_info(pset)
  search_fields <- pgx_candidate_columns(
    info,
    c("sample", "cell", "name", "alias", "id", "tissue", "histology", "source")
  )
  matched_rows <- pgx_find_rows(info, samples, match_mode, search_fields)
  result <- pgx_select_fields(
    info[matched_rows, , drop = FALSE],
    "sample",
    fields
  )
  result$matched_fields <- pgx_matched_fields(
    info,
    matched_rows,
    samples,
    match_mode,
    search_fields
  )
  limit <- max(1L, as.integer(limit))

  list(
    dataset = dataset,
    query = pgx_clean_vector(samples),
    match_mode = pgx_resolve_match_mode(match_mode),
    total_matches = nrow(result),
    returned = min(nrow(result), limit),
    metadata = pgx_limit_records(result, limit),
    searched_fields = unique(c("sample", search_fields))
  )
}

pgx_get_treatment_metadata <- function(
  dataset = "GDSCsmall",
  treatments,
  fields = NULL,
  match_mode = "exact",
  limit = 25
) {
  pset <- pgx_load_dataset(dataset)
  info <- pgx_treatment_info(pset)
  search_fields <- pgx_candidate_columns(
    info,
    c(
      "treatment",
      "drug",
      "compound",
      "name",
      "synonym",
      "target",
      "mechanism",
      "class",
      "pubchem",
      "cid",
      "inchikey",
      "chembl"
    )
  )
  match_mode <- pgx_resolve_match_mode(match_mode)
  resolution <- if (identical(match_mode, "exact")) {
    pgx_resolve_annotation_queries(info, treatments, search_fields)
  } else {
    NULL
  }
  matched_rows <- if (is.null(resolution)) {
    pgx_find_rows(info, treatments, match_mode, search_fields)
  } else {
    resolution$rows
  }
  result <- pgx_select_fields(
    info[matched_rows, , drop = FALSE],
    "treatment",
    fields
  )
  result$matched_fields <- if (is.null(resolution)) {
    pgx_matched_fields(
      info,
      matched_rows,
      treatments,
      match_mode,
      search_fields
    )
  } else {
    rep(paste(resolution$matched_fields, collapse = ";"), nrow(result))
  }
  limit <- max(1L, as.integer(limit))

  list(
    dataset = dataset,
    query = pgx_clean_vector(treatments),
    match_mode = match_mode,
    total_matches = nrow(result),
    returned = min(nrow(result), limit),
    metadata = pgx_limit_records(result, limit),
    searched_fields = unique(c("treatment", search_fields))
  )
}

pgx_filter_samples <- function(
  dataset = "GDSCsmall",
  filters,
  operator = "and",
  match_mode = "contains",
  fields = NULL,
  limit = 50
) {
  pset <- pgx_load_dataset(dataset)
  info <- pgx_sample_info(pset)
  filters <- pgx_parse_filters(filters)
  operator <- match.arg(operator, c("and", "or"))
  match_mode <- pgx_resolve_match_mode(match_mode)
  limit <- max(1L, as.integer(limit))

  missing_fields <- setdiff(names(filters), colnames(info))
  if (length(missing_fields) > 0) {
    stop(
      "Unknown sample filter field(s): ",
      paste(missing_fields, collapse = ", "),
      ". Use pgx_list_available_covariates() to discover sample fields."
    )
  }

  filter_matrix <- vapply(
    names(filters),
    function(field) {
      pgx_match_filter_values(info[[field]], filters[[field]], match_mode)
    },
    logical(nrow(info))
  )
  if (is.null(dim(filter_matrix))) {
    filter_matrix <- matrix(filter_matrix, ncol = 1)
  }

  keep <- if (identical(operator, "and")) {
    rowSums(filter_matrix) == ncol(filter_matrix)
  } else {
    rowSums(filter_matrix) > 0
  }

  result <- pgx_select_fields(info[keep, , drop = FALSE], "sample", fields)
  result$matched_filters <- apply(
    filter_matrix[keep, , drop = FALSE],
    1,
    function(matches) paste(names(filters)[matches], collapse = ";")
  )

  list(
    dataset = dataset,
    filters = filters,
    operator = operator,
    match_mode = match_mode,
    total_matches = nrow(result),
    returned = min(nrow(result), limit),
    samples = pgx_limit_records(result, limit)
  )
}

pgx_compare_metrics <- function(
  dataset = "GDSCsmall",
  drug = NULL,
  metrics = NULL,
  top_n = 10,
  rank_direction = "auto"
) {
  pset <- pgx_load_dataset(dataset)
  available <- pgx_sensitivity_measures(pset)
  if (length(available) == 0) {
    stop(
      "No classic sensitivity measures were available. For combo PSets, use",
      " pgx_list_treatment_response_tables() and combo-specific tools."
    )
  }
  metrics <- pgx_clean_vector(metrics)
  if (is.null(metrics)) {
    metrics <- grep("auc|ic50", available, value = TRUE, ignore.case = TRUE)
    if (length(metrics) < 2) {
      metrics <- available
    }
  }

  missing_metrics <- setdiff(metrics, available)
  if (length(missing_metrics) > 0) {
    stop(
      "Unknown sensitivity measure(s): ",
      paste(missing_metrics, collapse = ", "),
      ". Available measures are: ",
      paste(available, collapse = ", ")
    )
  }
  if (length(metrics) < 2) {
    stop("At least two sensitivity measures are required for comparison.")
  }

  top_n <- max(1L, as.integer(top_n))
  rank_direction <- match.arg(rank_direction, c("auto", "lowest", "highest"))
  warnings <- character()
  records_by_metric <- lapply(metrics, function(metric) {
    pgx_sensitivity_records(pset, metric, drug)
  })
  names(records_by_metric) <- metrics

  joined <- Reduce(
    function(x, y) merge(x, y, by = c("treatment", "sample"), all = FALSE),
    records_by_metric
  )
  finite_counts <- do.call(
    rbind,
    lapply(metrics, function(metric) {
      values <- joined[[metric]]
      data.frame(
        metric = metric,
        total = length(values),
        missing = sum(is.na(values)),
        finite = sum(!is.na(values) & is.finite(values)),
        non_finite = sum(!is.na(values) & !is.finite(values)),
        stringsAsFactors = FALSE
      )
    })
  )

  if (top_n >= nrow(joined)) {
    warnings <- c(
      warnings,
      paste(
        "top_n is greater than or equal to the number of complete records;",
        "top-N overlap is not meaningful."
      )
    )
  }

  pairwise <- do.call(
    rbind,
    lapply(utils::combn(metrics, 2, simplify = FALSE), function(pair) {
      pair_values <- joined[, pair, drop = FALSE]
      complete <- stats::complete.cases(pair_values) &
        rowSums(!is.finite(as.matrix(pair_values))) == 0
      direction_1 <- pgx_metric_direction(pair[[1]], rank_direction)
      direction_2 <- pgx_metric_direction(pair[[2]], rank_direction)
      expected_sign <- if (identical(direction_1, direction_2)) {
        "positive"
      } else {
        "negative"
      }
      spearman <- if (sum(complete) >= 2) {
        stats::cor(
          joined[[pair[[1]]]][complete],
          joined[[pair[[2]]]][complete],
          method = "spearman"
        )
      } else {
        NA_real_
      }
      unexpected <- !is.na(spearman) &&
        ((identical(expected_sign, "positive") && spearman < 0) ||
          (identical(expected_sign, "negative") && spearman > 0))
      data.frame(
        metric_1 = pair[[1]],
        metric_2 = pair[[2]],
        complete_pairs = sum(complete),
        expected_sign = expected_sign,
        spearman = spearman,
        pearson = if (sum(complete) >= 2) {
          stats::cor(
            joined[[pair[[1]]]][complete],
            joined[[pair[[2]]]][complete],
            method = "pearson"
          )
        } else {
          NA_real_
        },
        unexpected_direction = unexpected,
        stringsAsFactors = FALSE
      )
    })
  )
  if (any(pairwise$unexpected_direction, na.rm = TRUE)) {
    warnings <- c(
      warnings,
      "One or more metric correlations have an unexpected direction."
    )
  }

  ranked_lists <- lapply(metrics, function(metric) {
    records <- records_by_metric[[metric]]
    direction <- pgx_metric_direction(metric, rank_direction)
    records$item_id <- if (pgx_is_blank(drug)) {
      paste(records$treatment, records$sample, sep = "::")
    } else {
      records$sample
    }
    sensitive <- records[
      order(
        records[[metric]],
        decreasing = identical(direction, "highest")
      ),
      ,
      drop = FALSE
    ]
    resistant <- records[
      order(
        records[[metric]],
        decreasing = !identical(direction, "highest")
      ),
      ,
      drop = FALSE
    ]

    list(
      metric = metric,
      direction = direction,
      sensitive = head(sensitive$item_id, top_n),
      resistant = head(resistant$item_id, top_n)
    )
  })
  names(ranked_lists) <- metrics

  overlap <- do.call(
    rbind,
    lapply(utils::combn(metrics, 2, simplify = FALSE), function(pair) {
      first <- ranked_lists[[pair[[1]]]]
      second <- ranked_lists[[pair[[2]]]]
      data.frame(
        metric_1 = pair[[1]],
        metric_2 = pair[[2]],
        top_sensitive_overlap = length(intersect(
          first$sensitive,
          second$sensitive
        )),
        top_resistant_overlap = length(intersect(
          first$resistant,
          second$resistant
        )),
        top_n = top_n,
        stringsAsFactors = FALSE
      )
    })
  )

  top_sensitive_records <- do.call(
    rbind,
    lapply(metrics, function(metric) {
      records <- records_by_metric[[metric]]
      direction <- pgx_metric_direction(metric, rank_direction)
      pgx_rank_metric_records(records, metric, direction, top_n, "value")
    })
  )
  top_resistant_records <- do.call(
    rbind,
    lapply(metrics, function(metric) {
      records <- records_by_metric[[metric]]
      direction <- pgx_metric_direction(metric, rank_direction)
      resistant_direction <- if (identical(direction, "highest")) {
        "lowest"
      } else {
        "highest"
      }
      pgx_rank_metric_records(
        records,
        metric,
        resistant_direction,
        top_n,
        "value"
      )
    })
  )

  list(
    dataset = dataset,
    drug = if (pgx_is_blank(drug)) NA_character_ else drug,
    metrics = metrics,
    compared_records = nrow(joined),
    rank_direction = rank_direction,
    metric_directions = data.frame(
      metric = metrics,
      direction = vapply(
        metrics,
        pgx_metric_direction,
        character(1),
        rank_direction
      ),
      stringsAsFactors = FALSE
    ),
    finite_counts = finite_counts,
    correlations = pairwise,
    overlap = overlap,
    top_sensitive_records = top_sensitive_records,
    top_resistant_records = top_resistant_records,
    warnings = unique(warnings),
    note = paste(
      "This is a descriptive robustness comparison. It is not a formal",
      "statistical validation."
    )
  )
}

pgx_get_dose_response_points <- function(
  dataset = "GDSCsmall",
  drug,
  sample,
  summarize_replicates = TRUE
) {
  pset <- pgx_load_dataset(dataset)
  response <- pgx_extract_dose_response(
    pset = pset,
    drug = drug,
    sample = sample,
    summarize_replicates = summarize_replicates
  )

  list(
    dataset = dataset,
    drug = drug,
    sample = sample,
    summarize_replicates = isTRUE(summarize_replicates),
    sensitivity_metadata = response$sensitivity_metadata,
    points = response$points,
    raw_points = response$raw_points
  )
}

pgx_plot_dose_response <- function(
  dataset = "GDSCsmall",
  drug,
  sample,
  output_dir = NULL,
  summarize_replicates = TRUE,
  fit_curve = TRUE
) {
  if (pgx_is_blank(output_dir)) {
    output_dir <- tempdir()
  }
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  response <- pgx_get_dose_response_points(
    dataset = dataset,
    drug = drug,
    sample = sample,
    summarize_replicates = summarize_replicates
  )
  points <- response$points
  if (!all(c("dose", "viability") %in% colnames(points))) {
    stop("Dose-response points must include dose and viability columns.")
  }

  metrics <- pgx_compute_dose_response_metrics(
    concentration = points$dose,
    viability = points$viability,
    viability_as_pct = TRUE,
    area_type = "Actual"
  )

  fit <- NULL
  if (isTRUE(fit_curve) && nrow(points) >= 3) {
    fit_input <- data.frame(
      cell_id = sample,
      drug_id = drug,
      conc = points$dose,
      viability = points$viability / 100
    )
    fit <- pgx_capture(PharmacoGx::curveFittingPGX(
      fit_input,
      output_type = "all",
      main_fit_func = "hill"
    ))
  }

  output_file <- file.path(
    output_dir,
    paste0(pgx_safe_filename("PharmacoGx", dataset, drug, sample), ".png")
  )

  grDevices::png(output_file, width = 1000, height = 800, res = 120)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::plot(
    points$dose,
    points$viability,
    log = "x",
    pch = 19,
    xlab = "Dose",
    ylab = "Viability (%)",
    main = paste(dataset, drug, sample, sep = " / "),
    ylim = range(c(0, 100, points$viability), na.rm = TRUE)
  )
  if (!is.null(fit) && isTRUE(fit$ok) && !is.null(fit$value$curves)) {
    curves <- fit$value$curves
    fitted_viability <- curves$fitted_viability
    if (max(fitted_viability, na.rm = TRUE) <= 1) {
      fitted_viability <- fitted_viability * 100
    }
    graphics::lines(
      curves$conc,
      fitted_viability,
      col = "firebrick",
      lwd = 2
    )
    graphics::legend(
      "topright",
      legend = c("Observed", "Hill fit"),
      pch = c(19, NA),
      lty = c(NA, 1),
      col = c("black", "firebrick"),
      bty = "n"
    )
  }

  list(
    dataset = dataset,
    drug = drug,
    sample = sample,
    output_file = normalizePath(output_file, mustWork = FALSE),
    points = points,
    metrics = metrics$metrics,
    fitted_metrics = if (!is.null(fit) && isTRUE(fit$ok)) {
      fit$value$metrics
    } else {
      data.frame()
    },
    warnings = if (!is.null(fit)) fit$warnings else character()
  )
}

pgx_pset_curation_questions <- function(workflow = "unknown") {
  workflow <- match.arg(
    workflow,
    c("unknown", "sensitivity", "combination", "perturbation")
  )

  common_questions <- data.frame(
    topic = c("samples", "treatments", "responses", "doses", "metadata"),
    question = c(
      "Which table identifies samples or cell lines, and what is the sample ID column?",
      "Which table identifies treatments or compounds, and what is the treatment ID column?",
      "Which table contains response or viability measurements?",
      "Which columns contain dose values and dose units?",
      "Which sample, treatment, tissue, lineage, batch, or assay metadata are available?"
    ),
    stringsAsFactors = FALSE
  )

  combo_questions <- data.frame(
    topic = c("combination_treatments", "dose_matrix", "monotherapy_controls"),
    question = c(
      "For combinations, which columns identify treatment 1 and treatment 2?",
      "Which columns identify treatment 1 dose and treatment 2 dose for each matrix point?",
      "Are matched monotherapy and untreated controls present in the same table?"
    ),
    stringsAsFactors = FALSE
  )

  list(
    workflow = workflow,
    questions = if (identical(workflow, "combination")) {
      rbind(common_questions, combo_questions)
    } else {
      common_questions
    },
    next_step = paste(
      "Call pgx_validate_pset_inputs() with the available column names before",
      "attempting PharmacoSet construction."
    )
  )
}

pgx_validate_pset_inputs <- function(
  workflow = "sensitivity",
  sample_columns = NULL,
  treatment_columns = NULL,
  response_columns = NULL,
  dose_columns = NULL,
  metadata_columns = NULL
) {
  workflow <- match.arg(
    workflow,
    c("sensitivity", "combination", "perturbation")
  )

  required <- switch(
    workflow,
    "sensitivity" = list(
      sample = c("sampleid"),
      treatment = c("treatmentid"),
      response = c("sampleid", "treatmentid", "viability"),
      dose = c("sampleid", "treatmentid", "dose")
    ),
    "combination" = list(
      sample = c("sampleid"),
      treatment = c("treatment1id", "treatment2id"),
      response = c("sampleid", "treatment1id", "treatment2id", "viability"),
      dose = c("treatment1dose", "treatment2dose")
    ),
    "perturbation" = list(
      sample = c("sampleid"),
      treatment = c("treatmentid"),
      response = c("sampleid", "treatmentid"),
      dose = c("dose")
    )
  )

  results <- rbind(
    pgx_required_column_result("sample", sample_columns, required$sample),
    pgx_required_column_result(
      "treatment",
      treatment_columns,
      required$treatment
    ),
    pgx_required_column_result("response", response_columns, required$response),
    pgx_required_column_result("dose", dose_columns, required$dose)
  )

  list(
    workflow = workflow,
    ok = all(results$ok),
    validation = results,
    metadata_columns = pgx_clean_vector(metadata_columns),
    recommendation = if (all(results$ok)) {
      "Required column names are present. Next step is row-level validation and ID consistency checks."
    } else {
      "Add or map the missing columns before attempting PharmacoSet construction."
    }
  )
}

pgx_get_molecular_profile <- function(
  dataset = "GDSCsmall",
  mDataType,
  features,
  samples = NULL,
  feature_match_mode = "exact",
  summary_stat = NULL,
  limit_features = 25,
  limit_samples = 25,
  limit = NULL
) {
  if (!pgx_is_blank(limit)) {
    limit_features <- as.integer(limit)
    limit_samples <- as.integer(limit)
  }
  pset <- pgx_load_dataset(dataset)
  available_profiles <- pgx_mdata_names(pset)
  if (!mDataType %in% available_profiles) {
    return(list(
      dataset = dataset,
      mDataType = mDataType,
      available = FALSE,
      available_molecular_profiles = available_profiles,
      matched_features = data.frame(),
      dimensions = c(features = 0L, samples = 0L),
      records = pgx_matrix_records(matrix(numeric(), nrow = 0, ncol = 0)),
      note = paste0(
        "Molecular profile '",
        mDataType,
        "' is not available. Call pgx_list_available_covariates() before ",
        "requesting molecular data."
      )
    ))
  }
  capabilities <- pgx_molecular_profile_capabilities(pset)
  selected_capability <- capabilities[
    capabilities$molecular_profile == mDataType,
    ,
    drop = FALSE
  ]
  if (
    nrow(selected_capability) == 1 &&
      !isTRUE(selected_capability$usable[[1]])
  ) {
    return(list(
      dataset = dataset,
      mDataType = mDataType,
      available = FALSE,
      available_molecular_profiles = available_profiles,
      molecular_profile_capability = selected_capability,
      matched_features = data.frame(),
      dimensions = c(features = 0L, samples = 0L),
      records = pgx_matrix_records(matrix(numeric(), nrow = 0, ncol = 0)),
      note = paste(
        "The molecular profile slot exists but contains no usable feature-by-",
        "sample matrix. Do not treat the profile name alone as data availability."
      )
    ))
  }

  molecular <- pgx_get_molecular_matrix(
    pset = pset,
    mDataType = mDataType,
    features = features,
    samples = samples,
    feature_match_mode = feature_match_mode,
    summary_stat = summary_stat,
    limit_features = limit_features,
    limit_samples = limit_samples
  )

  list(
    dataset = dataset,
    mDataType = mDataType,
    available = TRUE,
    available_molecular_profiles = available_profiles,
    summary_stat = molecular$summary_stat,
    matched_features = molecular$matched_features,
    dimensions = stats::setNames(
      as.integer(dim(molecular$matrix)),
      c("features", "samples")
    ),
    records = pgx_matrix_records(molecular$matrix),
    note = if (nrow(molecular$matched_features) == 0) {
      "No requested features matched this molecular profile."
    } else {
      "Values are bounded to the requested feature and sample limits."
    }
  )
}

pgx_association_test <- function(
  dataset = "GDSCsmall",
  drug,
  metric = "auc_recomputed",
  mDataType,
  features,
  method = "spearman",
  samples = NULL,
  feature_match_mode = "exact",
  summary_stat = NULL,
  p_adjust_method = "BH"
) {
  method <- match.arg(method, c("spearman", "pearson", "wilcoxon"))
  pset <- pgx_load_dataset(dataset)
  available_metrics <- pgx_sensitivity_measures(pset)
  metric_resolution <- pgx_resolve_metric(metric, available_metrics)
  if (is.null(metric_resolution$resolved)) {
    return(list(
      dataset = dataset,
      drug = drug,
      metric = metric,
      mDataType = mDataType,
      available = FALSE,
      available_sensitivity_measures = available_metrics,
      suggested_sensitivity_measures = metric_resolution$suggestions,
      matched_features = data.frame(),
      results = data.frame(),
      note = paste0(
        "Sensitivity measure '",
        metric,
        "' is unavailable. Inspect available_sensitivity_measures and retry."
      )
    ))
  }
  metric <- metric_resolution$resolved
  response <- pgx_summarize_sensitivity_profiles(
    pset = pset,
    sensitivity.measure = metric,
    drugs = drug
  )
  treatment_rows <- pgx_summary_treatment_rows(pset, response, drug)
  if (length(treatment_rows) != 1) {
    return(list(
      dataset = dataset,
      drug = drug,
      metric = metric,
      mDataType = mDataType,
      available = FALSE,
      available_sensitivity_measures = available_metrics,
      resolved_treatments = treatment_rows,
      matched_features = data.frame(),
      results = data.frame(),
      note = if (length(treatment_rows) == 0) {
        paste0(
          "Treatment '",
          drug,
          "' did not match the stored treatment identifiers."
        )
      } else {
        paste0(
          "Treatment '",
          drug,
          "' matched multiple stored identifiers. Retry with one exact ID."
        )
      }
    ))
  }

  response_values <- response[treatment_rows[[1]], , drop = TRUE]
  if (!is.null(samples)) {
    response_values <- response_values[names(response_values) %in% samples]
  }

  molecular <- pgx_get_molecular_matrix(
    pset = pset,
    mDataType = mDataType,
    features = features,
    samples = names(response_values),
    feature_match_mode = feature_match_mode,
    summary_stat = summary_stat,
    limit_features = 100,
    limit_samples = length(response_values)
  )
  mat <- molecular$matrix
  if (nrow(mat) == 0 || ncol(mat) == 0) {
    return(list(
      dataset = dataset,
      drug = drug,
      metric = metric,
      mDataType = mDataType,
      matched_features = molecular$matched_features,
      results = data.frame(),
      note = paste(
        "No matching molecular features were found. Exploratory association",
        "only; bundled small datasets are demo fixtures."
      )
    ))
  }

  results <- pgx_feature_response_association(
    mat = mat,
    response_values = response_values,
    method = method,
    p_adjust_method = p_adjust_method
  )

  list(
    dataset = dataset,
    drug = drug,
    resolved_treatment = treatment_rows[[1]],
    metric = metric,
    mDataType = mDataType,
    available = TRUE,
    matched_features = molecular$matched_features,
    results = results,
    note = paste(
      "Exploratory association only. Bundled small datasets are demo fixtures",
      "and are not biologically powered for biomarker claims."
    )
  )
}

pgx_suggest_biomarker_candidates <- function(drug = NULL, limit = 20) {
  candidates <- pgx_known_biomarker_candidates
  candidates$candidate_source <- "curated workflow hypothesis"
  candidates$availability_checked <- FALSE
  drug <- pgx_clean_vector(drug)
  limit <- max(1L, as.integer(limit))

  if (!is.null(drug)) {
    searchable <- paste(candidates$drug, candidates$drug_aliases, sep = ";")
    keep <- pgx_match_filter_values(searchable, drug, "contains")
    candidates <- candidates[keep, , drop = FALSE]
  }

  list(
    query = drug,
    total = nrow(candidates),
    returned = min(nrow(candidates), limit),
    candidates = pgx_limit_records(candidates, limit),
    note = paste(
      "These are workflow seed candidates, not evidence and not confirmation",
      "that a feature exists in a PSet. Call pgx_check_analysis_requirements()",
      "before pgx_association_test() or pgx_multi_pset_association_test()."
    )
  )
}

pgx_multi_pset_association_test <- function(
  datasets = c("GDSCsmall", "CCLEsmall"),
  drug,
  metric = "auc_recomputed",
  mDataType,
  features,
  method = "spearman",
  feature_match_mode = "exact",
  summary_stat = NULL,
  p_adjust_method = "BH"
) {
  datasets <- pgx_clean_vector(datasets)
  if (is.null(datasets) || length(datasets) < 2) {
    stop("At least two datasets or local PSet file paths are required.")
  }
  method <- match.arg(method, c("spearman", "pearson", "wilcoxon"))

  per_dataset <- list()
  errors <- list()
  for (dataset in datasets) {
    result <- pgx_capture(pgx_association_test(
      dataset = dataset,
      drug = drug,
      metric = metric,
      mDataType = mDataType,
      features = features,
      method = method,
      feature_match_mode = feature_match_mode,
      summary_stat = summary_stat,
      p_adjust_method = p_adjust_method
    ))

    if (isTRUE(result$ok) && !isFALSE(result$value$available)) {
      records <- result$value$results
      if (nrow(records) > 0) {
        records$dataset <- dataset
        per_dataset[[dataset]] <- records[,
          c(
            "dataset",
            setdiff(colnames(records), "dataset")
          ),
          drop = FALSE
        ]
      }
    } else {
      error_message <- if (isTRUE(result$ok)) {
        result$value$note
      } else {
        result$value$message
      }
      errors[[dataset]] <- data.frame(
        dataset = dataset,
        error = error_message,
        stringsAsFactors = FALSE
      )
    }
  }

  per_dataset <- if (length(per_dataset) > 0) {
    do.call(rbind, per_dataset)
  } else {
    data.frame()
  }
  errors <- if (length(errors) > 0) {
    do.call(rbind, errors)
  } else {
    data.frame()
  }

  consensus <- if (nrow(per_dataset) > 0) {
    valid <- per_dataset[
      !is.na(per_dataset$p_value) &
        is.finite(per_dataset$p_value) &
        !is.na(per_dataset$effect_size) &
        is.finite(per_dataset$effect_size),
      ,
      drop = FALSE
    ]

    if (nrow(valid) > 0) {
      do.call(
        rbind,
        lapply(split(valid, valid$feature), function(feature_rows) {
          directions <- ifelse(
            feature_rows$effect_size > 0,
            "positive",
            "negative"
          )
          positive <- sum(directions == "positive")
          negative <- sum(directions == "negative")
          supported <- max(positive, negative)
          data.frame(
            feature = feature_rows$feature[[1]],
            datasets_tested = length(unique(per_dataset$dataset[
              per_dataset$feature == feature_rows$feature[[1]]
            ])),
            datasets_with_p_value = nrow(feature_rows),
            positive_effects = positive,
            negative_effects = negative,
            consensus_direction = if (positive >= negative) {
              "positive"
            } else {
              "negative"
            },
            direction_support = supported / max(1L, positive + negative),
            median_effect_size = stats::median(
              feature_rows$effect_size,
              na.rm = TRUE
            ),
            min_p_value = min(feature_rows$p_value, na.rm = TRUE),
            fisher_p_value = pgx_fisher_p(feature_rows$p_value),
            stringsAsFactors = FALSE
          )
        })
      )
    } else {
      data.frame()
    }
  } else {
    data.frame()
  }

  if (nrow(consensus) > 0) {
    consensus$fisher_fdr <- stats::p.adjust(
      consensus$fisher_p_value,
      method = p_adjust_method
    )
    consensus <- consensus[
      order(consensus$fisher_fdr, -consensus$direction_support),
      ,
      drop = FALSE
    ]
  }

  list(
    datasets = datasets,
    drug = drug,
    metric = metric,
    mDataType = mDataType,
    features = pgx_clean_vector(features),
    per_dataset_results = per_dataset,
    consensus = consensus,
    errors = errors,
    note = paste(
      "Each PSet was analyzed independently. Raw response or molecular data",
      "were not pooled across PSets because of batch-effect risk."
    )
  )
}

pgx_list_treatment_response_tables <- function(
  dataset = "GDSCsmall",
  limit_columns = 50
) {
  pset <- pgx_load_dataset(dataset)
  treatment_response <- pgx_treatment_response(pset)
  limit_columns <- max(1L, as.integer(limit_columns))

  tables <- do.call(
    rbind,
    lapply(names(treatment_response), function(table_name) {
      x <- treatment_response[[table_name]]
      columns <- colnames(x)
      data.frame(
        table = table_name,
        rows = nrow(x),
        columns = ncol(x),
        id_columns = paste(
          intersect(
            c(
              "treatmentid",
              "treatment1id",
              "treatment2id",
              "sampleid",
              "cellid",
              "tech_rep",
              "plateid",
              "BARCODE"
            ),
            columns
          ),
          collapse = ", "
        ),
        score_candidates = paste(pgx_score_candidates(x), collapse = ", "),
        returned_columns = paste(
          pgx_trim(columns, limit_columns),
          collapse = ", "
        ),
        stringsAsFactors = FALSE
      )
    })
  )

  list(
    dataset = dataset,
    tables = tables,
    note = paste(
      "Use score_candidates with combo ranking or combo biomarker tools. Local",
      "combo PSets often store useful endpoints in treatmentResponse profiles",
      "or synergy tables rather than classic sensitivity matrices."
    )
  )
}

pgx_rank_synergy_combinations <- function(
  dataset,
  table = "auto",
  score = NULL,
  rank_direction = "auto",
  min_samples = 3,
  limit = 20
) {
  .datatable.aware <- TRUE
  pset <- pgx_load_dataset(dataset)
  response <- pgx_response_table(
    pset,
    table = table,
    preferred = c("combo_viability", "synergy", "profiles", "raw"),
    as_data_table = TRUE
  )
  x <- data.table::as.data.table(response$table)
  score <- pgx_choose_score_column(x, score)
  rank_direction <- pgx_score_direction(score, rank_direction)
  min_samples <- max(1L, as.integer(min_samples))
  limit <- max(1L, as.integer(limit))

  required <- c("treatment1id", "treatment2id", "sampleid")
  missing <- setdiff(required, colnames(x))
  if (length(missing) > 0) {
    stop(
      "Combo ranking requires column(s): ",
      paste(missing, collapse = ", ")
    )
  }

  x <- x[
    !is.na(x[[score]]) & is.finite(x[[score]]),
    ,
    drop = FALSE
  ]
  sample_scores <- x[,
    list(
      sample_score = stats::median(get(score), na.rm = TRUE),
      records = .N
    ),
    by = c("treatment1id", "treatment2id", "sampleid")
  ]
  combo_scores <- sample_scores[,
    list(
      n_samples = data.table::uniqueN(sampleid),
      n_records = sum(records),
      mean_score = mean(sample_score, na.rm = TRUE),
      median_score = stats::median(sample_score, na.rm = TRUE),
      min_score = min(sample_score, na.rm = TRUE),
      max_score = max(sample_score, na.rm = TRUE)
    ),
    by = c("treatment1id", "treatment2id")
  ]
  combo_scores <- combo_scores[n_samples >= min_samples]
  data.table::setorderv(
    combo_scores,
    "median_score",
    order = if (identical(rank_direction, "highest")) -1L else 1L
  )
  combo_scores <- combo_scores[seq_len(min(nrow(combo_scores), limit))]
  combo_scores$rank <- seq_len(nrow(combo_scores))
  data.table::setcolorder(combo_scores, "rank")

  list(
    dataset = dataset,
    table = response$name,
    score = score,
    rank_direction = rank_direction,
    min_samples = min_samples,
    total_ranked = nrow(combo_scores),
    combinations = as.data.frame(combo_scores),
    note = paste(
      "Scores are aggregated per sample first, then ranked across combinations.",
      "Direction is inferred from the score name unless rank_direction is set."
    )
  )
}

pgx_get_combo_response_records <- function(
  dataset,
  treatment1 = NULL,
  treatment2 = NULL,
  table = "auto",
  score = NULL,
  match_mode = "exact",
  symmetric = TRUE,
  fields = NULL,
  limit = 100
) {
  pset <- pgx_load_dataset(dataset)
  response <- pgx_response_table(
    pset,
    table = table,
    preferred = c("combo_viability", "synergy", "profiles", "raw"),
    as_data_table = FALSE
  )
  x <- pgx_filter_combo_rows(
    response$table,
    treatment1 = treatment1,
    treatment2 = treatment2,
    match_mode = match_mode,
    symmetric = symmetric
  )
  limit <- max(1L, as.integer(limit))

  score <- if (pgx_is_blank(score)) {
    candidates <- pgx_score_candidates(x)
    if (length(candidates) > 0) candidates[[1]] else NULL
  } else {
    pgx_choose_score_column(x, score)
  }

  default_fields <- unique(c(
    "treatment1id",
    "treatment2id",
    "treatment1dose",
    "treatment2dose",
    "sampleid",
    "tech_rep",
    "plateid",
    "BARCODE",
    score,
    pgx_trim(pgx_score_candidates(x), 10)
  ))
  fields <- pgx_clean_vector(fields)
  missing_fields <- character()
  if (is.null(fields)) {
    keep_fields <- intersect(default_fields, colnames(x))
  } else {
    missing_fields <- setdiff(fields, colnames(x))
    keep_fields <- intersect(fields, colnames(x))
    if (length(keep_fields) == 0) {
      keep_fields <- intersect(default_fields, colnames(x))
    }
  }

  list(
    dataset = dataset,
    table = response$name,
    treatment1 = pgx_clean_vector(treatment1),
    treatment2 = pgx_clean_vector(treatment2),
    match_mode = pgx_resolve_match_mode(match_mode),
    symmetric = isTRUE(symmetric),
    score = score,
    total_matches = nrow(x),
    returned = min(nrow(x), limit),
    requested_fields = fields,
    unavailable_fields = missing_fields,
    available_fields = colnames(x),
    note = if (length(missing_fields) > 0) {
      paste0(
        "Unavailable requested fields were omitted: ",
        paste(missing_fields, collapse = ", "),
        ". Inspect available_fields before retrying."
      )
    } else {
      "Returned fields are present in the selected combination-response table."
    },
    records = pgx_limit_records(x[, keep_fields, drop = FALSE], limit)
  )
}

pgx_combo_biomarker_association <- function(
  dataset,
  treatment1,
  treatment2,
  mDataType,
  features,
  table = "auto",
  score = NULL,
  method = "spearman",
  match_mode = "exact",
  symmetric = TRUE,
  feature_match_mode = "exact",
  summary_stat = NULL,
  p_adjust_method = "BH"
) {
  .datatable.aware <- TRUE
  method <- match.arg(method, c("spearman", "pearson", "wilcoxon"))
  pset <- pgx_load_dataset(dataset)
  response <- pgx_response_table(
    pset,
    table = table,
    preferred = c("combo_viability", "synergy", "profiles", "raw"),
    as_data_table = TRUE
  )
  x <- pgx_filter_combo_rows(
    response$table,
    treatment1 = treatment1,
    treatment2 = treatment2,
    match_mode = match_mode,
    symmetric = symmetric
  )
  x <- data.table::as.data.table(x)
  score <- pgx_choose_score_column(x, score)

  if (!"sampleid" %in% colnames(x)) {
    stop("Combo biomarker association requires a sampleid column.")
  }
  x <- x[
    !is.na(x[[score]]) & is.finite(x[[score]]),
    ,
    drop = FALSE
  ]
  if (nrow(x) == 0) {
    stop("No finite combo phenotype values remained after filtering.")
  }

  phenotype <- x[,
    list(
      phenotype = stats::median(get(score), na.rm = TRUE),
      records = .N
    ),
    by = "sampleid"
  ]
  response_values <- phenotype$phenotype
  names(response_values) <- phenotype$sampleid

  molecular <- pgx_get_molecular_matrix(
    pset = pset,
    mDataType = mDataType,
    features = features,
    samples = names(response_values),
    feature_match_mode = feature_match_mode,
    summary_stat = summary_stat,
    limit_features = 100,
    limit_samples = length(response_values)
  )
  results <- pgx_feature_response_association(
    mat = molecular$matrix,
    response_values = response_values,
    method = method,
    p_adjust_method = p_adjust_method
  )

  list(
    dataset = dataset,
    table = response$name,
    treatment1 = pgx_clean_vector(treatment1),
    treatment2 = pgx_clean_vector(treatment2),
    score = score,
    method = method,
    mDataType = mDataType,
    matched_features = molecular$matched_features,
    phenotype_summary = data.frame(
      samples = length(response_values),
      records = sum(phenotype$records),
      min = min(response_values, na.rm = TRUE),
      median = stats::median(response_values, na.rm = TRUE),
      max = max(response_values, na.rm = TRUE),
      stringsAsFactors = FALSE
    ),
    results = results,
    note = paste(
      "Combo phenotype values are median-aggregated per sample before",
      "feature-response association. Interpret as exploratory biomarker",
      "support for this combo endpoint."
    )
  )
}

pgx_summarize_numeric <- function(x, summary_stat = "median") {
  x <- x[is.finite(x)]
  if (length(x) == 0) {
    return(NA_real_)
  }

  switch(
    summary_stat,
    median = stats::median(x, na.rm = TRUE),
    mean = mean(x, na.rm = TRUE),
    max = max(x, na.rm = TRUE),
    min = min(x, na.rm = TRUE),
    stop("Unsupported summary_stat: ", summary_stat)
  )
}

pgx_mono_biomarker_association_for_pset <- function(
  pset,
  treatment,
  mDataType,
  features,
  mono_table = "mono_profiles",
  mono_metric = "aac_recomputed",
  method = "spearman",
  match_mode = "exact",
  feature_match_mode = "exact",
  summary_stat = NULL,
  p_adjust_method = "BH"
) {
  .datatable.aware <- TRUE
  response <- pgx_response_table(
    pset,
    table = mono_table,
    preferred = c("mono_profiles", "profiles", "raw"),
    as_data_table = TRUE
  )
  x <- response$table
  x <- data.table::as.data.table(x)
  treatment_col <- intersect(c("treatmentid", "treatment1id"), colnames(x))
  if (length(treatment_col) == 0 || !"sampleid" %in% colnames(x)) {
    stop(
      "Mono biomarker comparison requires treatmentid or treatment1id plus sampleid."
    )
  }
  treatment_col <- treatment_col[[1]]
  if (!mono_metric %in% colnames(x) || !is.numeric(x[[mono_metric]])) {
    stop(
      "Mono metric '",
      mono_metric,
      "' is not a numeric column in treatmentResponse table '",
      response$name,
      "'."
    )
  }

  keep <- pgx_match_filter_values(x[[treatment_col]], treatment, match_mode)
  x <- x[keep & is.finite(x[[mono_metric]]), , drop = FALSE]
  if (nrow(x) == 0) {
    return(list(
      table = response$name,
      treatment = treatment,
      metric = mono_metric,
      phenotype_summary = data.frame(),
      matched_features = data.frame(),
      results = data.frame()
    ))
  }

  phenotype <- x[,
    list(
      phenotype = stats::median(get(mono_metric), na.rm = TRUE),
      records = .N
    ),
    by = "sampleid"
  ]
  response_values <- phenotype$phenotype
  names(response_values) <- phenotype$sampleid

  molecular <- pgx_get_molecular_matrix(
    pset = pset,
    mDataType = mDataType,
    features = features,
    samples = names(response_values),
    feature_match_mode = feature_match_mode,
    summary_stat = summary_stat,
    limit_features = 100,
    limit_samples = length(response_values)
  )
  results <- pgx_feature_response_association(
    mat = molecular$matrix,
    response_values = response_values,
    method = method,
    p_adjust_method = p_adjust_method
  )

  list(
    table = response$name,
    treatment = treatment,
    metric = mono_metric,
    phenotype_summary = data.frame(
      samples = length(response_values),
      records = sum(phenotype$records),
      min = min(response_values, na.rm = TRUE),
      median = stats::median(response_values, na.rm = TRUE),
      max = max(response_values, na.rm = TRUE),
      stringsAsFactors = FALSE
    ),
    matched_features = molecular$matched_features,
    results = results
  )
}

pgx_prefix_result_columns <- function(results, prefix) {
  if (nrow(results) == 0) {
    return(data.frame(feature = character(), stringsAsFactors = FALSE))
  }
  out <- results
  names(out)[names(out) != "feature"] <- paste0(
    prefix,
    "_",
    names(out)[names(out) != "feature"]
  )
  out
}

pgx_effect_direction <- function(x) {
  direction <- rep(NA_character_, length(x))
  direction[is.finite(x) & x > 0] <- "positive"
  direction[is.finite(x) & x < 0] <- "negative"
  direction[is.finite(x) & x == 0] <- "zero"
  direction
}

pgx_direction_agreement <- function(x, y) {
  out <- rep(NA_character_, length(x))
  comparable <- !is.na(x) & !is.na(y) & x != "zero" & y != "zero"
  out[comparable] <- ifelse(x[comparable] == y[comparable], "same", "opposite")
  out
}

pgx_compare_mono_combo_biomarkers <- function(
  dataset,
  treatment1,
  treatment2,
  mDataType,
  features,
  combo_table = "auto",
  combo_score = NULL,
  mono_table = "mono_profiles",
  mono_metric = "aac_recomputed",
  method = "spearman",
  match_mode = "exact",
  symmetric = TRUE,
  feature_match_mode = "exact",
  summary_stat = NULL,
  p_adjust_method = "BH"
) {
  .datatable.aware <- TRUE
  method <- match.arg(method, c("spearman", "pearson", "wilcoxon"))
  pset <- pgx_load_dataset(dataset)

  combo_response <- pgx_response_table(
    pset,
    table = combo_table,
    preferred = c("combo_viability", "synergy", "profiles", "raw"),
    as_data_table = TRUE
  )
  combo_rows <- pgx_filter_combo_rows(
    combo_response$table,
    treatment1 = treatment1,
    treatment2 = treatment2,
    match_mode = match_mode,
    symmetric = symmetric
  )
  combo_rows <- data.table::as.data.table(combo_rows)
  combo_score <- pgx_choose_score_column(combo_rows, combo_score)
  combo_rows <- combo_rows[
    is.finite(combo_rows[[combo_score]]),
    ,
    drop = FALSE
  ]
  if (!"sampleid" %in% colnames(combo_rows) || nrow(combo_rows) == 0) {
    stop("No finite combo phenotype rows with sample IDs remained.")
  }

  combo_phenotype <- combo_rows[,
    list(
      phenotype = stats::median(get(combo_score), na.rm = TRUE),
      records = .N
    ),
    by = "sampleid"
  ]
  combo_response_values <- combo_phenotype$phenotype
  names(combo_response_values) <- combo_phenotype$sampleid
  combo_molecular <- pgx_get_molecular_matrix(
    pset = pset,
    mDataType = mDataType,
    features = features,
    samples = names(combo_response_values),
    feature_match_mode = feature_match_mode,
    summary_stat = summary_stat,
    limit_features = 100,
    limit_samples = length(combo_response_values)
  )
  combo_results <- pgx_feature_response_association(
    mat = combo_molecular$matrix,
    response_values = combo_response_values,
    method = method,
    p_adjust_method = p_adjust_method
  )

  mono1 <- pgx_mono_biomarker_association_for_pset(
    pset = pset,
    treatment = treatment1,
    mDataType = mDataType,
    features = features,
    mono_table = mono_table,
    mono_metric = mono_metric,
    method = method,
    match_mode = match_mode,
    feature_match_mode = feature_match_mode,
    summary_stat = summary_stat,
    p_adjust_method = p_adjust_method
  )
  mono2 <- pgx_mono_biomarker_association_for_pset(
    pset = pset,
    treatment = treatment2,
    mDataType = mDataType,
    features = features,
    mono_table = mono_table,
    mono_metric = mono_metric,
    method = method,
    match_mode = match_mode,
    feature_match_mode = feature_match_mode,
    summary_stat = summary_stat,
    p_adjust_method = p_adjust_method
  )

  comparison <- Reduce(
    function(x, y) merge(x, y, by = "feature", all = TRUE),
    list(
      pgx_prefix_result_columns(combo_results, "combo"),
      pgx_prefix_result_columns(mono1$results, "mono1"),
      pgx_prefix_result_columns(mono2$results, "mono2")
    )
  )
  comparison$combo_direction <- pgx_effect_direction(
    comparison$combo_effect_size
  )
  comparison$mono1_direction <- pgx_effect_direction(
    comparison$mono1_effect_size
  )
  comparison$mono2_direction <- pgx_effect_direction(
    comparison$mono2_effect_size
  )
  comparison$combo_vs_mono1_direction <- pgx_direction_agreement(
    comparison$combo_direction,
    comparison$mono1_direction
  )
  comparison$combo_vs_mono2_direction <- pgx_direction_agreement(
    comparison$combo_direction,
    comparison$mono2_direction
  )

  list(
    dataset = dataset,
    treatment1 = treatment1,
    treatment2 = treatment2,
    mDataType = mDataType,
    features = pgx_clean_vector(features),
    method = method,
    combo_table = combo_response$name,
    combo_score = combo_score,
    mono_table = mono1$table,
    mono_metric = mono_metric,
    matched_features = combo_molecular$matched_features,
    phenotype_summaries = list(
      combo = data.frame(
        samples = length(combo_response_values),
        records = sum(combo_phenotype$records),
        min = min(combo_response_values, na.rm = TRUE),
        median = stats::median(combo_response_values, na.rm = TRUE),
        max = max(combo_response_values, na.rm = TRUE),
        stringsAsFactors = FALSE
      ),
      mono1 = mono1$phenotype_summary,
      mono2 = mono2$phenotype_summary
    ),
    comparison = comparison,
    note = paste(
      "This bounded workflow compares feature-response associations for one",
      "combo phenotype against each monotherapy arm. It is exploratory and",
      "does not run genome-wide biomarker discovery unless the caller supplies",
      "a broad feature list."
    )
  )
}

pgx_plot_synergy_heatmap <- function(
  dataset,
  treatment1 = NULL,
  treatment2 = NULL,
  sample = NULL,
  table = "auto",
  score = NULL,
  rank_direction = "auto",
  match_mode = "exact",
  symmetric = TRUE,
  output_dir = NULL,
  summary_stat = "median"
) {
  .datatable.aware <- TRUE
  summary_stat <- match.arg(summary_stat, c("median", "mean", "max", "min"))
  if (pgx_is_blank(output_dir)) {
    output_dir <- tempdir()
  }
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  pset <- pgx_load_dataset(dataset)
  response <- pgx_response_table(
    pset,
    table = table,
    preferred = c("combo_viability", "synergy", "raw", "profiles"),
    as_data_table = TRUE
  )
  x <- data.table::as.data.table(response$table)
  score <- pgx_choose_score_column(x, score)
  rank_direction <- pgx_score_direction(score, rank_direction)

  if (pgx_is_blank(treatment1) || pgx_is_blank(treatment2)) {
    required <- c("treatment1id", "treatment2id", "sampleid")
    missing <- setdiff(required, colnames(x))
    if (length(missing) > 0) {
      stop(
        "Auto-selecting a combo requires column(s): ",
        paste(missing, collapse = ", ")
      )
    }
    ranked <- x[is.finite(x[[score]]), ][,
      list(
        sample_score = stats::median(get(score), na.rm = TRUE),
        records = .N
      ),
      by = c("treatment1id", "treatment2id", "sampleid")
    ][,
      list(
        n_samples = data.table::uniqueN(sampleid),
        median_score = stats::median(sample_score, na.rm = TRUE),
        n_records = sum(records)
      ),
      by = c("treatment1id", "treatment2id")
    ]
    data.table::setorderv(
      ranked,
      "median_score",
      order = if (identical(rank_direction, "highest")) -1L else 1L
    )
    if (nrow(ranked) == 0) {
      stop("No finite combo scores were available for heatmap selection.")
    }
    treatment1 <- ranked$treatment1id[[1]]
    treatment2 <- ranked$treatment2id[[1]]
  }

  x <- pgx_filter_combo_rows(
    x,
    treatment1 = treatment1,
    treatment2 = treatment2,
    match_mode = match_mode,
    symmetric = symmetric
  )
  x <- data.table::as.data.table(x)
  if (!pgx_is_blank(sample)) {
    if (!"sampleid" %in% colnames(x)) {
      stop("Sample filtering requires a sampleid column.")
    }
    x <- x[
      pgx_match_filter_values(x$sampleid, sample, match_mode),
      ,
      drop = FALSE
    ]
  }
  if (nrow(x) == 0) {
    stop("No combo rows matched the requested treatments/sample.")
  }

  required <- c("treatment1dose", "treatment2dose", "sampleid")
  missing <- setdiff(required, colnames(x))
  if (length(missing) > 0) {
    stop(
      "Synergy heatmap requires column(s): ",
      paste(missing, collapse = ", ")
    )
  }
  x <- x[
    is.finite(x[[score]]) &
      is.finite(x[["treatment1dose"]]) &
      is.finite(x[["treatment2dose"]]),
    ,
    drop = FALSE
  ]
  if (nrow(x) == 0) {
    stop("No finite dose/score rows remained for the heatmap.")
  }

  if (pgx_is_blank(sample)) {
    sample_scores <- x[,
      list(sample_score = stats::median(get(score), na.rm = TRUE)),
      by = "sampleid"
    ]
    data.table::setorderv(
      sample_scores,
      "sample_score",
      order = if (identical(rank_direction, "highest")) -1L else 1L
    )
    sample <- sample_scores$sampleid[[1]]
    x <- x[x[["sampleid"]] == sample, , drop = FALSE]
  }

  grid <- x[,
    list(
      score_value = pgx_summarize_numeric(get(score), summary_stat),
      records = .N
    ),
    by = c("treatment1dose", "treatment2dose")
  ]
  dose1 <- sort(unique(grid$treatment1dose))
  dose2 <- sort(unique(grid$treatment2dose))
  z <- matrix(
    NA_real_,
    nrow = length(dose1),
    ncol = length(dose2),
    dimnames = list(dose1, dose2)
  )
  for (i in seq_len(nrow(grid))) {
    z[
      match(grid$treatment1dose[[i]], dose1),
      match(grid$treatment2dose[[i]], dose2)
    ] <- grid$score_value[[i]]
  }

  output_file <- file.path(
    output_dir,
    paste0(
      pgx_safe_filename(
        "PharmacoGx_synergy_heatmap",
        dataset,
        treatment1,
        treatment2,
        sample,
        score
      ),
      ".png"
    )
  )

  z_range <- range(z, na.rm = TRUE)
  if (!all(is.finite(z_range))) {
    stop("Heatmap score matrix contains no finite values.")
  }
  z_abs <- max(abs(z_range), na.rm = TRUE)
  colors <- grDevices::colorRampPalette(c("#2c7bb6", "white", "#d7191c"))(101)

  grDevices::png(output_file, width = 1100, height = 900, res = 120)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::image(
    x = seq_along(dose1),
    y = seq_along(dose2),
    z = z,
    col = colors,
    zlim = c(-z_abs, z_abs),
    axes = FALSE,
    xlab = "Treatment 1 dose",
    ylab = "Treatment 2 dose",
    main = paste(
      dataset,
      treatment1,
      treatment2,
      sample,
      score,
      sep = " / "
    )
  )
  graphics::axis(
    1,
    at = seq_along(dose1),
    labels = signif(dose1, 3),
    las = 2,
    cex.axis = 0.75
  )
  graphics::axis(
    2,
    at = seq_along(dose2),
    labels = signif(dose2, 3),
    las = 2,
    cex.axis = 0.75
  )
  graphics::box()

  list(
    dataset = dataset,
    table = response$name,
    treatment1 = as.character(treatment1),
    treatment2 = as.character(treatment2),
    sample = as.character(sample),
    score = score,
    rank_direction = rank_direction,
    summary_stat = summary_stat,
    output_file = normalizePath(output_file, mustWork = FALSE),
    grid = as.data.frame(grid),
    note = paste(
      "Heatmap cells summarize matched combo rows at each dose pair. Positive",
      "values indicate synergy for standard expected-minus-observed score",
      "columns such as Bliss_score, HSA_score, ZIP_score, or ZIP_delta."
    )
  )
}

pgx_detect_assay_mode <- function(dataset = "GDSCsmall") {
  pset <- pgx_load_dataset(dataset)
  info <- pgx_sensitivity_info(pset)
  treatment_names <- tryCatch(
    pgx_treatment_names(pset),
    error = function(e) character()
  )
  treatment_response_tables <- tryCatch(
    pgx_treatment_response(pset),
    error = function(e) {
      list()
    }
  )
  dataset_type <- tryCatch(PharmacoGx::datasetType(pset), error = function(e) {
    character()
  })
  has_sensitivity <- nrow(info) > 0
  has_combo_cols <- any(c("treatment2id", "treatment2dose") %in% colnames(info))
  combo_rows <- if ("treatment2id" %in% colnames(info)) {
    sum(!is.na(info$treatment2id) & nzchar(as.character(info$treatment2id)))
  } else {
    0L
  }
  combo_names <- sum(grepl("///", treatment_names, fixed = TRUE))
  combo_table_rows <- sum(vapply(
    names(treatment_response_tables),
    function(table_name) {
      x <- treatment_response_tables[[table_name]]
      all(c("treatment1id", "treatment2id") %in% colnames(x)) && nrow(x) > 0
    },
    logical(1)
  ))
  has_combo <- has_combo_cols &&
    combo_rows > 0 ||
    combo_names > 0 ||
    combo_table_rows > 0
  has_mono <- has_sensitivity && (nrow(info) > combo_rows)

  mode <- if (!has_sensitivity && "perturbation" %in% dataset_type) {
    "perturbation_only"
  } else if (has_combo && has_mono) {
    "mixed"
  } else if (has_combo) {
    "combination"
  } else if (has_sensitivity) {
    "monotherapy"
  } else {
    "unknown"
  }

  list(
    dataset = dataset,
    dataset_type = dataset_type,
    assay_mode = mode,
    sensitivity_rows = nrow(info),
    combo_rows = combo_rows,
    combo_treatment_names = combo_names,
    combo_treatment_response_tables = combo_table_rows,
    note = if (has_combo) {
      "Combination-like fields were detected; inspect raw data before synergy analysis."
    } else {
      "No combination-dose fields were detected in this dataset."
    }
  )
}

pgx_find_pset_overlap <- function(
  datasets = c("GDSCsmall", "CCLEsmall"),
  limit = 50
) {
  psets <- pgx_load_many_datasets(datasets)
  limit <- max(1L, as.integer(limit))

  sample_sets <- lapply(psets, pgx_sample_names)
  treatment_sets <- lapply(psets, pgx_treatment_names)
  pair_sets <- lapply(psets, pgx_pair_table)

  shared_samples <- Reduce(intersect, sample_sets)
  shared_treatments <- Reduce(intersect, treatment_sets)
  pair_keys <- lapply(pair_sets, function(x) {
    paste(x$sample, x$treatment, sep = "::")
  })
  shared_pair_keys <- Reduce(intersect, pair_keys)
  shared_pairs <- if (length(shared_pair_keys) > 0) {
    do.call(
      rbind,
      strsplit(shared_pair_keys, "::", fixed = TRUE)
    )
  } else {
    matrix(character(), ncol = 2)
  }
  shared_pairs <- data.frame(
    sample = shared_pairs[, 1],
    treatment = shared_pairs[, 2],
    stringsAsFactors = FALSE
  )

  list(
    datasets = names(psets),
    sample_counts = vapply(sample_sets, length, integer(1)),
    treatment_counts = vapply(treatment_sets, length, integer(1)),
    pair_counts = vapply(pair_sets, nrow, integer(1)),
    shared_sample_count = length(shared_samples),
    shared_treatment_count = length(shared_treatments),
    shared_pair_count = nrow(shared_pairs),
    shared_samples = pgx_trim(shared_samples, limit),
    shared_treatments = pgx_trim(shared_treatments, limit),
    shared_pairs = pgx_limit_records(shared_pairs, limit)
  )
}

pgx_compare_pset_response <- function(
  datasets = c("GDSCsmall", "CCLEsmall"),
  drug,
  metric = "auc_recomputed",
  samples = NULL,
  limit = 50
) {
  psets <- pgx_load_many_datasets(datasets)
  samples <- pgx_clean_vector(samples)
  limit <- max(1L, as.integer(limit))

  records <- do.call(
    rbind,
    lapply(names(psets), function(dataset) {
      if (!metric %in% pgx_sensitivity_measures(psets[[dataset]])) {
        return(data.frame())
      }
      mat <- pgx_summarize_sensitivity_profiles(
        pset = psets[[dataset]],
        sensitivity.measure = metric,
        drugs = drug
      )
      treatment_rows <- pgx_summary_treatment_rows(
        psets[[dataset]],
        mat,
        drug
      )
      if (length(treatment_rows) != 1) {
        return(data.frame())
      }
      values <- mat[treatment_rows[[1]], , drop = TRUE]
      if (!is.null(samples)) {
        values <- values[names(values) %in% samples]
      }
      data.frame(
        dataset = dataset,
        sample = names(values),
        requested_treatment = drug,
        treatment = treatment_rows[[1]],
        metric = metric,
        value = unname(as.numeric(values)),
        stringsAsFactors = FALSE
      )
    })
  )

  if (nrow(records) == 0) {
    return(list(
      datasets = names(psets),
      drug = drug,
      metric = metric,
      records = data.frame(),
      correlations = data.frame()
    ))
  }

  wide <- reshape(
    records[, c("dataset", "sample", "value")],
    idvar = "sample",
    timevar = "dataset",
    direction = "wide"
  )
  value_cols <- grep("^value\\.", colnames(wide), value = TRUE)
  correlations <- if (length(value_cols) >= 2) {
    do.call(
      rbind,
      lapply(utils::combn(value_cols, 2, simplify = FALSE), function(pair) {
        complete <- stats::complete.cases(wide[, pair, drop = FALSE])
        data.frame(
          dataset_1 = sub("^value\\.", "", pair[[1]]),
          dataset_2 = sub("^value\\.", "", pair[[2]]),
          complete_samples = sum(complete),
          spearman = if (sum(complete) >= 2) {
            stats::cor(
              wide[[pair[[1]]]][complete],
              wide[[pair[[2]]]][complete],
              method = "spearman"
            )
          } else {
            NA_real_
          },
          stringsAsFactors = FALSE
        )
      })
    )
  } else {
    data.frame()
  }

  list(
    datasets = names(psets),
    drug = drug,
    metric = metric,
    total_records = nrow(records),
    records = pgx_limit_records(records, limit),
    correlations = correlations,
    note = "Identifier matching is exact; inspect metadata before biological interpretation."
  )
}

pgx_summarize_sensitivity <- function(
  dataset = "GDSCsmall",
  sensitivity_measure = "auc_recomputed",
  drugs = NULL,
  cell_lines = NULL,
  limit = 50
) {
  pset <- pgx_load_dataset(dataset)
  drugs <- pgx_clean_vector(drugs)
  cell_lines <- pgx_clean_vector(cell_lines)
  limit <- max(1L, as.integer(limit))
  available_metrics <- pgx_sensitivity_measures(pset)
  if (!sensitivity_measure %in% available_metrics) {
    return(list(
      dataset = dataset,
      sensitivity_measure = sensitivity_measure,
      available = FALSE,
      available_sensitivity_measures = available_metrics,
      available_response_tables = pgx_treatment_response_table_names(pset),
      dimensions = c(treatments = 0L, samples = 0L),
      records = data.frame(),
      note = paste0(
        "Sensitivity measure '",
        sensitivity_measure,
        "' is unavailable. Inspect available_sensitivity_measures and retry."
      )
    ))
  }

  summary <- pgx_summarize_sensitivity_profiles(
    pset = pset,
    sensitivity.measure = sensitivity_measure,
    drugs = drugs,
    cell.lines = cell_lines
  )

  list(
    dataset = dataset,
    sensitivity_measure = sensitivity_measure,
    available = TRUE,
    requested_treatments = drugs,
    resolved_treatments = attr(
      summary,
      "resolved_treatments",
      exact = TRUE
    ),
    response_table = attr(summary, "response_table", exact = TRUE),
    dimensions = stats::setNames(
      as.integer(dim(summary)),
      c("treatments", "samples")
    ),
    records = pgx_matrix_to_records(
      summary,
      value_name = "sensitivity",
      limit = limit
    )
  )
}

pgx_top_responders <- function(
  dataset = "GDSCsmall",
  drug,
  sensitivity_measure = "auc_recomputed",
  rank_direction = "auto",
  limit = 10
) {
  rank_direction <- match.arg(rank_direction, c("auto", "lowest", "highest"))
  limit <- max(1L, as.integer(limit))
  pset <- pgx_load_dataset(dataset)
  available_metrics <- pgx_sensitivity_measures(pset)
  if (!sensitivity_measure %in% available_metrics) {
    return(list(
      dataset = dataset,
      drug = drug,
      sensitivity_measure = sensitivity_measure,
      available = FALSE,
      available_sensitivity_measures = available_metrics,
      top_responders = data.frame(),
      note = paste0(
        "Sensitivity measure '",
        sensitivity_measure,
        "' is unavailable. Inspect available_sensitivity_measures and retry."
      )
    ))
  }

  summary <- pgx_summarize_sensitivity_profiles(
    pset = pset,
    sensitivity.measure = sensitivity_measure,
    drugs = drug
  )

  treatment_rows <- pgx_summary_treatment_rows(pset, summary, drug)
  if (length(treatment_rows) != 1) {
    return(list(
      dataset = dataset,
      drug = drug,
      sensitivity_measure = sensitivity_measure,
      available = FALSE,
      available_sensitivity_measures = available_metrics,
      resolved_treatments = treatment_rows,
      top_responders = data.frame(),
      note = if (length(treatment_rows) == 0) {
        paste0(
          "Treatment '",
          drug,
          "' did not match the stored treatment identifiers."
        )
      } else {
        paste0(
          "Treatment '",
          drug,
          "' matched multiple stored identifiers. Retry with one exact ID."
        )
      }
    ))
  }

  values <- summary[treatment_rows[[1]], , drop = TRUE]
  values <- values[!is.na(values)]

  if (identical(rank_direction, "auto")) {
    rank_direction <- if (
      grepl("auc|ic50", sensitivity_measure, ignore.case = TRUE)
    ) {
      "lowest"
    } else {
      "highest"
    }
  }

  order_index <- order(
    values,
    decreasing = identical(rank_direction, "highest")
  )
  values <- values[order_index]
  values <- values[seq_len(min(length(values), limit))]

  result <- data.frame(
    rank = seq_along(values),
    sample = names(values),
    value = unname(as.numeric(values)),
    stringsAsFactors = FALSE
  )

  list(
    dataset = dataset,
    drug = drug,
    resolved_treatment = treatment_rows[[1]],
    sensitivity_measure = sensitivity_measure,
    available = TRUE,
    rank_direction = rank_direction,
    interpretation = paste(
      "For AUC and IC50 measures, lower values usually indicate stronger",
      "sensitivity. For AAC-like response-area measures, higher values usually",
      "indicate stronger sensitivity."
    ),
    top_responders = result
  )
}

pgx_compute_dose_response_metrics <- function(
  concentration,
  viability,
  viability_as_pct = TRUE,
  area_type = "Actual"
) {
  if (length(concentration) != length(viability)) {
    stop("concentration and viability must have the same length.")
  }
  if (length(concentration) < 2) {
    stop("At least two concentration and viability values are required.")
  }

  area_type <- match.arg(area_type, c("Actual", "Fitted"))
  concentration <- as.numeric(concentration)
  viability <- as.numeric(viability)
  viability_as_pct <- isTRUE(viability_as_pct)

  list(
    input = data.frame(
      concentration = concentration,
      viability = viability
    ),
    metrics = list(
      auc = pgx_metric_value(PharmacoGx::computeAUC(
        concentration,
        viability,
        viability_as_pct = viability_as_pct,
        area.type = area_type,
        verbose = FALSE
      )),
      aac = pgx_metric_value(PharmacoGx::computeAAC(
        concentration,
        viability,
        viability_as_pct = viability_as_pct,
        area.type = area_type,
        verbose = FALSE
      )),
      ic50 = pgx_metric_value(PharmacoGx::computeIC50(
        concentration,
        viability,
        viability_as_pct = viability_as_pct,
        verbose = FALSE
      )),
      ac50 = pgx_metric_value(PharmacoGx::computeAC50(
        concentration,
        viability,
        viability_as_pct = viability_as_pct,
        verbose = FALSE
      ))
    ),
    note = paste(
      "AUC is normalized viability area in PharmacoGx 3.14.0 and later;",
      "AAC is the complementary response area."
    )
  )
}

pgx_compute_synergy_reference <- function(
  viability_1,
  viability_2,
  hsa_na_rm = FALSE
) {
  if (
    !length(viability_1) %in% c(1L, length(viability_2)) &&
      !length(viability_2) %in% c(1L, length(viability_1))
  ) {
    stop("viability_1 and viability_2 must have length 1 or the same length.")
  }

  viability_1 <- as.numeric(viability_1)
  viability_2 <- as.numeric(viability_2)
  hsa_arguments <- list(viability_1, viability_2)
  if ("na.rm" %in% names(formals(PharmacoGx::computeHSA))) {
    hsa_arguments$na.rm <- isTRUE(hsa_na_rm)
    hsa_reference <- do.call(PharmacoGx::computeHSA, hsa_arguments)
  } else {
    hsa_reference <- pmin(
      viability_1,
      viability_2,
      na.rm = isTRUE(hsa_na_rm)
    )
  }

  list(
    hsa_na_rm = isTRUE(hsa_na_rm),
    references = data.frame(
      viability_1 = viability_1,
      viability_2 = viability_2,
      bliss_reference = PharmacoGx::computeBliss(viability_1, viability_2),
      hsa_reference = hsa_reference
    )
  )
}

pgx_list_available_psets <- function(canonical = TRUE, limit = 25) {
  limit <- max(1L, as.integer(limit))
  available <- PharmacoGx::availablePSets(canonical = isTRUE(canonical))

  list(
    canonical = isTRUE(canonical),
    total = nrow(available),
    returned = min(nrow(available), limit),
    psets = pgx_limit_records(available, limit),
    note = paste(
      "Use the PSet Name value with pgx_download_pset(). Downloads may be",
      "large and require confirm_download = TRUE."
    )
  )
}

pgx_pset_summary <- function(pset, file_path = NA_character_) {
  capture_value <- function(expr, default = NA) {
    value <- tryCatch(expr, error = function(e) default)
    if (length(value) == 0) {
      return(default)
    }
    value
  }

  sensitivity_measures <- capture_value(
    pgx_sensitivity_measures(pset),
    character()
  )
  molecular_profiles <- capture_value(pgx_mdata_names(pset), character())

  list(
    name = capture_value(PharmacoGx::name(pset), NA_character_),
    class = paste(class(pset), collapse = ", "),
    file_path = file_path,
    dataset_type = capture_value(PharmacoGx::datasetType(pset), NA_character_),
    sample_count = capture_value(
      length(pgx_sample_names(pset)),
      NA_integer_
    ),
    treatment_count = capture_value(
      length(pgx_treatment_names(pset)),
      NA_integer_
    ),
    sensitivity_measures = sensitivity_measures,
    molecular_profiles = molecular_profiles,
    note = paste(
      "The downloaded PharmacoSet was not returned directly because S4 objects",
      "are not suitable MCP payloads. Use its saved file path in R for deeper",
      "analysis."
    )
  )
}

pgx_download_pset <- function(
  name,
  save_dir = NULL,
  pset_file_name = NULL,
  timeout = 600,
  confirm_download = FALSE
) {
  if (!isTRUE(confirm_download)) {
    return(list(
      downloaded = FALSE,
      confirmation_required = TRUE,
      requested_pset = name,
      message = paste(
        "Downloading a PharmacoSet can transfer large external files. Call this",
        "tool again with confirm_download = TRUE after the user explicitly",
        "approves the dataset name, runtime, storage location, and network use."
      )
    ))
  }

  if (pgx_is_blank(save_dir)) {
    save_dir <- file.path(
      tools::R_user_dir("PharmacoGx", which = "data"),
      "psets"
    )
  }
  dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

  if (pgx_is_blank(pset_file_name)) {
    pset_file_name <- NULL
  }

  file_path <- if (is.null(pset_file_name)) {
    file.path(save_dir, paste0(name, ".rds"))
  } else {
    file.path(save_dir, pset_file_name)
  }

  pset <- PharmacoGx::downloadPSet(
    name = name,
    saveDir = save_dir,
    pSetFileName = pset_file_name,
    verbose = FALSE,
    timeout = as.numeric(timeout)
  )
  if (!file.exists(file_path) || file.info(file_path)$size <= 0) {
    stop("The download call completed without a readable non-empty PSet file.")
  }

  list(
    downloaded = TRUE,
    requested_pset = name,
    save_dir = normalizePath(save_dir, mustWork = FALSE),
    file_path = normalizePath(file_path, mustWork = FALSE),
    local_alias = pgx_local_alias_from_path(file_path),
    automatically_discoverable = normalizePath(file_path, mustWork = TRUE) %in%
      pgx_discover_local_psets()$path_normalized,
    pset_summary = pgx_pset_summary(
      pset,
      file_path = normalizePath(file_path, mustWork = FALSE)
    )
  )
}

pgx_session_info <- function() {
  packages <- c("PharmacoGx", "CoreGx", "mcptools", "ellmer", "btw")
  versions <- vapply(
    packages,
    function(pkg) {
      if (!requireNamespace(pkg, quietly = TRUE)) {
        return(NA_character_)
      }
      as.character(utils::packageVersion(pkg))
    },
    character(1)
  )

  tool_file <- Sys.getenv("PHARMACOGX_MCP_TOOL_FILE", "")
  if (!nzchar(tool_file)) {
    tool_file <- normalizePath(
      file.path(Sys.getenv("PHARMACOGX_MCP_DIR", ""), "PharmacoGx-tools.R"),
      mustWork = FALSE
    )
  }
  tool_modified_at <- if (file.exists(tool_file)) {
    format(file.info(tool_file)$mtime, "%Y-%m-%dT%H:%M:%S%z")
  } else {
    NA_character_
  }
  server_started_at <- Sys.getenv(
    "PHARMACOGX_MCP_SERVER_STARTED_AT",
    NA_character_
  )
  tool_stale <- if (
    file.exists(tool_file) &&
      !is.na(server_started_at) &&
      nzchar(server_started_at)
  ) {
    file.info(tool_file)$mtime >
      as.POSIXct(
        server_started_at,
        format = "%Y-%m-%dT%H:%M:%S%z"
      )
  } else {
    NA
  }
  local_catalog <- tryCatch(
    pgx_local_pset_catalog(include_missing = FALSE),
    error = function(e) pgx_empty_local_manifest()
  )

  list(
    r_version = as.character(getRversion()),
    package_versions = data.frame(
      package = names(versions),
      version = unname(versions),
      installed = !is.na(versions),
      stringsAsFactors = FALSE
    ),
    working_directory = getwd(),
    supported_demo_datasets = pgx_supported_datasets,
    mcp_process_id = Sys.getpid(),
    mcp_server_started_at = server_started_at,
    pharmacogx_package_source = Sys.getenv(
      "PHARMACOGX_MCP_PACKAGE_SOURCE",
      "installed"
    ),
    mcp_tool_file = tool_file,
    mcp_tool_file_modified_at = tool_modified_at,
    restart_required = tool_stale,
    single_content_response = pgx_env_flag(
      "PHARMACOGX_MCP_SINGLE_CONTENT",
      default = TRUE
    ),
    mcp_log_file = pgx_log_file(),
    tool_count = if (exists("res", inherits = TRUE)) {
      length(res)
    } else {
      NA_integer_
    },
    local_manifest = normalizePath(
      pgx_default_local_manifest_path(),
      mustWork = FALSE
    ),
    local_discovery_roots = pgx_default_local_pset_roots(),
    discoverable_local_psets = local_catalog$name,
    cached_datasets = names(pgx_dataset_cache),
    capability_cache_directory = pgx_capability_cache_dir(),
    cached_capability_snapshots = names(pgx_capability_cache)
  )
}

if (!requireNamespace("ellmer", quietly = TRUE)) {
  stop(
    "The PharmacoGx MCP demo requires the ellmer package. Install it with ",
    "pak::pkg_install('ellmer') before starting mcptools::mcp_server()."
  )
}

if (!requireNamespace("PharmacoGx", quietly = TRUE)) {
  stop(
    "The PharmacoGx package must be installed before starting this MCP server."
  )
}

res <- list(
  pgx_list_example_datasets = pgx_tool(
    fun = pgx_list_example_datasets,
    name = "pgx_list_example_datasets",
    description = paste(
      "List bundled PharmacoGx demo datasets and the type of agentic",
      "analysis each supports."
    ),
    arguments = list()
  ),
  pgx_list_local_psets = pgx_tool(
    fun = pgx_list_local_psets,
    name = "pgx_list_local_psets",
    description = paste(
      "List local PSet manifest aliases such as nci_almanac, gdsc2_matrix,",
      "and gdsc2_anchor for real PSet-backed demos."
    ),
    arguments = list(
      manifest_path = ellmer::type_string(
        "Optional manifest CSV path. Defaults to inst/mcp/local-pset-manifest.csv.",
        required = FALSE
      ),
      roots = ellmer::type_array(
        "Optional directories to scan recursively for local .rds and .qs PSets.",
        items = ellmer::type_string(),
        required = FALSE
      ),
      discover = ellmer::type_boolean(
        "Whether to scan configured local PSet roots. Defaults to TRUE.",
        required = FALSE
      ),
      include_missing = ellmer::type_boolean(
        "Whether to include aliases whose files are missing. Defaults to TRUE.",
        required = FALSE
      )
    )
  ),
  pgx_check_analysis_requirements = pgx_tool(
    fun = pgx_check_analysis_requirements,
    name = "pgx_check_analysis_requirements",
    description = paste(
      "Check local and bundled PSets against required treatment, sample",
      "context, response metric, molecular profile, and feature. Resolves",
      "stored annotation identifiers such as PubChem CIDs and Ensembl IDs."
    ),
    arguments = list(
      datasets = ellmer::type_array(
        "Optional dataset names, aliases, or local paths to check.",
        items = ellmer::type_string(),
        required = FALSE
      ),
      drug = ellmer::type_string(
        "Optional treatment name or annotated identifier such as PubChem CID.",
        required = FALSE
      ),
      sample_query = ellmer::type_string(
        "Optional tissue, site, histology, lineage, disease, or panel query.",
        required = FALSE
      ),
      metric = ellmer::type_string(
        "Optional required sensitivity metric.",
        required = FALSE
      ),
      mDataType = ellmer::type_string(
        "Optional required molecular profile such as rna or rnaseq.",
        required = FALSE
      ),
      features = ellmer::type_array(
        "Optional required gene symbols, Ensembl IDs, or stored feature IDs.",
        items = ellmer::type_string(),
        required = FALSE
      ),
      include_local = ellmer::type_boolean(
        "Include auto-discovered local PSets when datasets are omitted.",
        required = FALSE
      ),
      include_bundled = ellmer::type_boolean(
        "Include bundled toy PSets when datasets are omitted.",
        required = FALSE
      ),
      limit = ellmer::type_integer(
        "Maximum number of datasets to inspect. Defaults to 20.",
        required = FALSE
      )
    )
  ),
  pgx_list_entities = pgx_tool(
    fun = pgx_list_entities,
    name = "pgx_list_entities",
    description = paste(
      "Search or page through treatments, samples, or sensitivity measures",
      "from a bundled or local PharmacoSet."
    ),
    arguments = list(
      dataset = ellmer::type_string(
        "One of GDSCsmall, CCLEsmall, or CMAPsmall. Defaults to GDSCsmall.",
        required = FALSE
      ),
      entity = ellmer::type_string(
        "Entity: treatments, samples, or sensitivity_measures.",
        required = FALSE
      ),
      query = ellmer::type_array(
        "Optional search values applied before pagination.",
        items = ellmer::type_string(),
        required = FALSE
      ),
      match_mode = ellmer::type_string(
        "One of exact or contains. Defaults to contains.",
        required = FALSE
      ),
      offset = ellmer::type_integer(
        "Zero-based result offset. Defaults to 0.",
        required = FALSE
      ),
      limit = ellmer::type_integer(
        "Maximum number of values to return. Defaults to 25.",
        required = FALSE
      )
    )
  ),
  pgx_list_available_covariates = pgx_tool(
    fun = pgx_list_available_covariates,
    name = "pgx_list_available_covariates",
    description = paste(
      "List sample, treatment, sensitivity, and molecular covariates",
      "available in a bundled or local PharmacoSet."
    ),
    arguments = list(
      dataset = ellmer::type_string(
        paste(
          "Bundled dataset, local manifest alias, or explicit .qs/.rds path.",
          "Defaults to GDSCsmall."
        ),
        required = FALSE
      )
    )
  ),
  pgx_get_sample_metadata = pgx_tool(
    fun = pgx_get_sample_metadata,
    name = "pgx_get_sample_metadata",
    description = paste(
      "Return sample or cell-line annotation rows from a bundled PharmacoGx",
      "dataset, including lineage and tissue fields when available."
    ),
    arguments = list(
      dataset = ellmer::type_string(
        "One of GDSCsmall, CCLEsmall, or CMAPsmall. Defaults to GDSCsmall.",
        required = FALSE
      ),
      samples = ellmer::type_array(
        "Sample or cell-line IDs/names to look up.",
        items = ellmer::type_string()
      ),
      fields = ellmer::type_array(
        "Optional metadata fields to return.",
        items = ellmer::type_string(),
        required = FALSE
      ),
      match_mode = ellmer::type_string(
        "One of exact or contains. Defaults to exact.",
        required = FALSE
      ),
      limit = ellmer::type_integer(
        "Maximum number of matching rows to return. Defaults to 25.",
        required = FALSE
      )
    )
  ),
  pgx_get_treatment_metadata = pgx_tool(
    fun = pgx_get_treatment_metadata,
    name = "pgx_get_treatment_metadata",
    description = paste(
      "Return treatment or drug annotation rows from a bundled PharmacoGx",
      "dataset, including cross-referenced identifiers, target, and mechanism",
      "fields when available."
    ),
    arguments = list(
      dataset = ellmer::type_string(
        "One of GDSCsmall, CCLEsmall, or CMAPsmall. Defaults to GDSCsmall.",
        required = FALSE
      ),
      treatments = ellmer::type_array(
        "Treatment names or identifiers such as PubChem CID or InChIKey.",
        items = ellmer::type_string()
      ),
      fields = ellmer::type_array(
        "Optional metadata fields to return.",
        items = ellmer::type_string(),
        required = FALSE
      ),
      match_mode = ellmer::type_string(
        "One of exact or contains. Defaults to exact.",
        required = FALSE
      ),
      limit = ellmer::type_integer(
        "Maximum number of matching rows to return. Defaults to 25.",
        required = FALSE
      )
    )
  ),
  pgx_filter_samples = pgx_tool(
    fun = pgx_filter_samples,
    name = "pgx_filter_samples",
    description = paste(
      "Filter sample or cell-line annotations by one or more fields from",
      "sampleInfo(). Call pgx_list_available_covariates first."
    ),
    arguments = list(
      dataset = ellmer::type_string(
        "One of GDSCsmall, CCLEsmall, or CMAPsmall. Defaults to GDSCsmall.",
        required = FALSE
      ),
      filters = ellmer::type_from_schema(
        text = paste0(
          '{"type":"object","description":"Named sample metadata filters, ',
          'for example {\\\"tissueid\\\": \\\"lung\\\"}.",',
          '"additionalProperties":true}'
        )
      ),
      operator = ellmer::type_string(
        "How to combine filters: and or or. Defaults to and.",
        required = FALSE
      ),
      match_mode = ellmer::type_string(
        "One of exact or contains. Defaults to contains.",
        required = FALSE
      ),
      fields = ellmer::type_array(
        "Optional metadata fields to return.",
        items = ellmer::type_string(),
        required = FALSE
      ),
      limit = ellmer::type_integer(
        "Maximum number of matching rows to return. Defaults to 50.",
        required = FALSE
      )
    )
  ),
  pgx_compare_metrics = pgx_tool(
    fun = pgx_compare_metrics,
    name = "pgx_compare_metrics",
    description = paste(
      "Compare sensitivity measures with rank correlations and top/bottom",
      "responder overlap for a drug or all treatment-sample pairs."
    ),
    arguments = list(
      dataset = ellmer::type_string(
        "One of GDSCsmall or CCLEsmall. Defaults to GDSCsmall.",
        required = FALSE
      ),
      drug = ellmer::type_string(
        "Optional treatment name. If omitted, compares all treatment-sample pairs.",
        required = FALSE
      ),
      metrics = ellmer::type_array(
        "Optional sensitivity measures to compare.",
        items = ellmer::type_string(),
        required = FALSE
      ),
      top_n = ellmer::type_integer(
        "Number of top sensitive/resistant records for overlap. Defaults to 10.",
        required = FALSE
      ),
      rank_direction = ellmer::type_string(
        "One of auto, lowest, or highest. Defaults to auto.",
        required = FALSE
      )
    )
  ),
  pgx_get_dose_response_points = pgx_tool(
    fun = pgx_get_dose_response_points,
    name = "pgx_get_dose_response_points",
    description = paste(
      "Retrieve raw concentration/viability points from sensitivityRaw()",
      "for a selected PharmacoGx sample-drug pair."
    ),
    arguments = list(
      dataset = ellmer::type_string(
        "One of GDSCsmall or CCLEsmall. Defaults to GDSCsmall.",
        required = FALSE
      ),
      drug = ellmer::type_string("Treatment or drug ID/name."),
      sample = ellmer::type_string("Sample or cell-line ID/name."),
      summarize_replicates = ellmer::type_boolean(
        "Whether to median summarize exact duplicate doses. Defaults to TRUE.",
        required = FALSE
      )
    )
  ),
  pgx_plot_dose_response = pgx_tool(
    fun = pgx_plot_dose_response,
    name = "pgx_plot_dose_response",
    description = paste(
      "Plot a PNG dose-response curve for a selected PharmacoGx sample-drug",
      "pair and return the file path plus metrics."
    ),
    arguments = list(
      dataset = ellmer::type_string(
        "One of GDSCsmall or CCLEsmall. Defaults to GDSCsmall.",
        required = FALSE
      ),
      drug = ellmer::type_string("Treatment or drug ID/name."),
      sample = ellmer::type_string("Sample or cell-line ID/name."),
      output_dir = ellmer::type_string(
        "Directory where the PNG should be written. Defaults to tempdir().",
        required = FALSE
      ),
      summarize_replicates = ellmer::type_boolean(
        "Whether to median summarize exact duplicate doses. Defaults to TRUE.",
        required = FALSE
      ),
      fit_curve = ellmer::type_boolean(
        "Whether to overlay a Hill fit when possible. Defaults to TRUE.",
        required = FALSE
      )
    )
  ),
  pgx_pset_curation_questions = pgx_tool(
    fun = pgx_pset_curation_questions,
    name = "pgx_pset_curation_questions",
    description = paste(
      "Return staged questions for collecting user table and column",
      "information before PharmacoSet curation."
    ),
    arguments = list(
      workflow = ellmer::type_string(
        "One of unknown, sensitivity, combination, or perturbation. Defaults to unknown.",
        required = FALSE
      )
    )
  ),
  pgx_validate_pset_inputs = pgx_tool(
    fun = pgx_validate_pset_inputs,
    name = "pgx_validate_pset_inputs",
    description = paste(
      "Validate whether supplied table column names cover the minimum",
      "PharmacoSet curation requirements for a workflow."
    ),
    arguments = list(
      workflow = ellmer::type_string(
        "One of sensitivity, combination, or perturbation. Defaults to sensitivity.",
        required = FALSE
      ),
      sample_columns = ellmer::type_array(
        "Columns present in the sample table.",
        items = ellmer::type_string(),
        required = FALSE
      ),
      treatment_columns = ellmer::type_array(
        "Columns present in the treatment table.",
        items = ellmer::type_string(),
        required = FALSE
      ),
      response_columns = ellmer::type_array(
        "Columns present in the response or viability table.",
        items = ellmer::type_string(),
        required = FALSE
      ),
      dose_columns = ellmer::type_array(
        "Columns present in the dose table or dose-response table.",
        items = ellmer::type_string(),
        required = FALSE
      ),
      metadata_columns = ellmer::type_array(
        "Optional metadata columns available for curation.",
        items = ellmer::type_string(),
        required = FALSE
      )
    )
  ),
  pgx_get_molecular_profile = pgx_tool(
    fun = pgx_get_molecular_profile,
    name = "pgx_get_molecular_profile",
    description = paste(
      "Return summarized molecular profile values for selected features and",
      "samples from a bundled or local PharmacoSet."
    ),
    arguments = list(
      dataset = ellmer::type_string(
        paste(
          "Bundled dataset, local manifest alias, or explicit .qs/.rds path.",
          "Defaults to GDSCsmall."
        ),
        required = FALSE
      ),
      mDataType = ellmer::type_string(
        "Molecular profile type such as rna, mutation, rnaseq, or cnv."
      ),
      features = ellmer::type_array(
        "Feature IDs or symbols to retrieve.",
        items = ellmer::type_string()
      ),
      samples = ellmer::type_array(
        "Optional sample IDs to retrieve.",
        items = ellmer::type_string(),
        required = FALSE
      ),
      feature_match_mode = ellmer::type_string(
        "One of exact or contains. Defaults to exact.",
        required = FALSE
      ),
      summary_stat = ellmer::type_string(
        "Optional summarizeMolecularProfiles summary.stat override.",
        required = FALSE
      ),
      limit_features = ellmer::type_integer(
        "Maximum number of matched features. Defaults to 25.",
        required = FALSE
      ),
      limit_samples = ellmer::type_integer(
        "Maximum number of samples. Defaults to 25.",
        required = FALSE
      ),
      limit = ellmer::type_integer(
        "Convenience limit applied to both matched features and samples.",
        required = FALSE
      )
    )
  ),
  pgx_association_test = pgx_tool(
    fun = pgx_association_test,
    name = "pgx_association_test",
    description = paste(
      "Run exploratory feature-response associations for one drug using",
      "summarized molecular profiles and sensitivity metrics.",
      "Spearman p-values are two-sided asymptotic (exact=FALSE).",
      "p_adjust_method applies within this call, not across separate datasets."
    ),
    arguments = list(
      dataset = ellmer::type_string(
        paste(
          "Bundled dataset, local manifest alias, or explicit .qs/.rds path.",
          "Defaults to GDSCsmall."
        ),
        required = FALSE
      ),
      drug = ellmer::type_string("Treatment or drug name."),
      metric = ellmer::type_string(
        "Sensitivity metric. Defaults to auc_recomputed.",
        required = FALSE
      ),
      mDataType = ellmer::type_string(
        "Molecular profile type such as rna or mutation."
      ),
      features = ellmer::type_array(
        "Feature IDs or symbols to test.",
        items = ellmer::type_string()
      ),
      method = ellmer::type_string(
        "One of spearman, pearson, or wilcoxon. Defaults to spearman.",
        required = FALSE
      ),
      samples = ellmer::type_array(
        "Optional sample IDs to include.",
        items = ellmer::type_string(),
        required = FALSE
      ),
      feature_match_mode = ellmer::type_string(
        "One of exact or contains. Defaults to exact.",
        required = FALSE
      ),
      summary_stat = ellmer::type_string(
        "Optional molecular summary.stat override.",
        required = FALSE
      ),
      p_adjust_method = ellmer::type_string(
        "p.adjust method. Defaults to BH.",
        required = FALSE
      )
    )
  ),
  pgx_suggest_biomarker_candidates = pgx_tool(
    fun = pgx_suggest_biomarker_candidates,
    name = "pgx_suggest_biomarker_candidates",
    description = paste(
      "Return curated seed biomarker candidates for agent workflows, such as",
      "SLC29A1/hENT1 for gemcitabine and PDE3A for anagrelide."
    ),
    arguments = list(
      drug = ellmer::type_string(
        "Optional drug name or alias to filter candidate biomarkers.",
        required = FALSE
      ),
      limit = ellmer::type_integer(
        "Maximum number of candidates to return. Defaults to 20.",
        required = FALSE
      )
    )
  ),
  pgx_multi_pset_association_test = pgx_tool(
    fun = pgx_multi_pset_association_test,
    name = "pgx_multi_pset_association_test",
    description = paste(
      "Run feature-response association independently across multiple",
      "datasets or local PSet file paths and summarize consensus evidence."
    ),
    arguments = list(
      datasets = ellmer::type_array(
        "Dataset names or explicit local .qs/.rds PSet file paths.",
        items = ellmer::type_string()
      ),
      drug = ellmer::type_string("Treatment or drug name."),
      metric = ellmer::type_string(
        "Sensitivity metric. Defaults to auc_recomputed.",
        required = FALSE
      ),
      mDataType = ellmer::type_string(
        "Molecular profile type such as rna, mutation, or cnv."
      ),
      features = ellmer::type_array(
        "Feature IDs or symbols to test.",
        items = ellmer::type_string()
      ),
      method = ellmer::type_string(
        "One of spearman, pearson, or wilcoxon. Defaults to spearman.",
        required = FALSE
      ),
      feature_match_mode = ellmer::type_string(
        "One of exact or contains. Defaults to exact.",
        required = FALSE
      ),
      summary_stat = ellmer::type_string(
        "Optional molecular summary.stat override.",
        required = FALSE
      ),
      p_adjust_method = ellmer::type_string(
        "p.adjust method. Defaults to BH.",
        required = FALSE
      )
    )
  ),
  pgx_list_treatment_response_tables = pgx_tool(
    fun = pgx_list_treatment_response_tables,
    name = "pgx_list_treatment_response_tables",
    description = paste(
      "List treatmentResponse tables, columns, and score candidates for a",
      "bundled or local PharmacoSet."
    ),
    arguments = list(
      dataset = ellmer::type_string(
        "Dataset name or explicit local .qs/.rds PSet file path."
      ),
      limit_columns = ellmer::type_integer(
        "Maximum number of column names to return per table. Defaults to 50.",
        required = FALSE
      )
    )
  ),
  pgx_rank_synergy_combinations = pgx_tool(
    fun = pgx_rank_synergy_combinations,
    name = "pgx_rank_synergy_combinations",
    description = paste(
      "Rank drug combinations from treatmentResponse combo tables using a",
      "selected or auto-detected synergy/response score."
    ),
    arguments = list(
      dataset = ellmer::type_string(
        "Dataset name or explicit local .qs/.rds combo PSet file path."
      ),
      table = ellmer::type_string(
        "treatmentResponse table name or auto. Defaults to auto.",
        required = FALSE
      ),
      score = ellmer::type_string(
        "Optional numeric score column. Auto-detected if omitted.",
        required = FALSE
      ),
      rank_direction = ellmer::type_string(
        "One of auto, lowest, or highest. Defaults to auto.",
        required = FALSE
      ),
      min_samples = ellmer::type_integer(
        "Minimum number of samples required per combo. Defaults to 3.",
        required = FALSE
      ),
      limit = ellmer::type_integer(
        "Maximum number of ranked combos to return. Defaults to 20.",
        required = FALSE
      )
    )
  ),
  pgx_get_combo_response_records = pgx_tool(
    fun = pgx_get_combo_response_records,
    name = "pgx_get_combo_response_records",
    description = paste(
      "Return bounded combo treatment-response rows for selected treatments",
      "from a bundled or local combo PharmacoSet."
    ),
    arguments = list(
      dataset = ellmer::type_string(
        "Dataset name or explicit local .qs/.rds combo PSet file path."
      ),
      treatment1 = ellmer::type_string(
        "Optional first treatment ID/name filter.",
        required = FALSE
      ),
      treatment2 = ellmer::type_string(
        "Optional second treatment ID/name filter.",
        required = FALSE
      ),
      table = ellmer::type_string(
        "treatmentResponse table name or auto. Defaults to auto.",
        required = FALSE
      ),
      score = ellmer::type_string(
        "Optional numeric score column to include.",
        required = FALSE
      ),
      match_mode = ellmer::type_string(
        "One of exact or contains. Defaults to exact.",
        required = FALSE
      ),
      symmetric = ellmer::type_boolean(
        "Whether treatment1/treatment2 order can be swapped. Defaults to TRUE.",
        required = FALSE
      ),
      fields = ellmer::type_array(
        paste(
          "Optional fields to return. Unknown fields are reported and omitted;",
          "inspect available_fields in the result."
        ),
        items = ellmer::type_string(),
        required = FALSE
      ),
      limit = ellmer::type_integer(
        "Maximum number of rows to return. Defaults to 100.",
        required = FALSE
      )
    )
  ),
  pgx_plot_synergy_heatmap = pgx_tool(
    fun = pgx_plot_synergy_heatmap,
    name = "pgx_plot_synergy_heatmap",
    description = paste(
      "Create a PNG heatmap for a selected or auto-selected combo dose matrix",
      "using a synergy or response score."
    ),
    arguments = list(
      dataset = ellmer::type_string(
        "Dataset alias, dataset name, or explicit local .qs/.rds combo PSet path."
      ),
      treatment1 = ellmer::type_string(
        "Optional first treatment ID/name. If omitted, the top combo is auto-selected.",
        required = FALSE
      ),
      treatment2 = ellmer::type_string(
        "Optional second treatment ID/name. If omitted, the top combo is auto-selected.",
        required = FALSE
      ),
      sample = ellmer::type_string(
        "Optional sample/cell-line ID. If omitted, the strongest sample is auto-selected.",
        required = FALSE
      ),
      table = ellmer::type_string(
        "treatmentResponse table name or auto. Defaults to auto.",
        required = FALSE
      ),
      score = ellmer::type_string(
        "Optional numeric score column such as ZIP_delta or Bliss_score.",
        required = FALSE
      ),
      rank_direction = ellmer::type_string(
        "One of auto, lowest, or highest. Defaults to auto.",
        required = FALSE
      ),
      match_mode = ellmer::type_string(
        "One of exact or contains. Defaults to exact.",
        required = FALSE
      ),
      symmetric = ellmer::type_boolean(
        "Whether treatment1/treatment2 order can be swapped. Defaults to TRUE.",
        required = FALSE
      ),
      output_dir = ellmer::type_string(
        "Directory where the PNG should be written. Defaults to tempdir().",
        required = FALSE
      ),
      summary_stat = ellmer::type_string(
        "One of median, mean, max, or min. Defaults to median.",
        required = FALSE
      )
    )
  ),
  pgx_combo_biomarker_association = pgx_tool(
    fun = pgx_combo_biomarker_association,
    name = "pgx_combo_biomarker_association",
    description = paste(
      "Associate selected molecular features with a combo response or synergy",
      "score after aggregating the combo phenotype per sample."
    ),
    arguments = list(
      dataset = ellmer::type_string(
        "Dataset name or explicit local .qs/.rds combo PSet file path."
      ),
      treatment1 = ellmer::type_string("First treatment ID/name."),
      treatment2 = ellmer::type_string("Second treatment ID/name."),
      mDataType = ellmer::type_string(
        "Molecular profile type such as mutation or cnv."
      ),
      features = ellmer::type_array(
        "Feature IDs or symbols to test.",
        items = ellmer::type_string()
      ),
      table = ellmer::type_string(
        "treatmentResponse table name or auto. Defaults to auto.",
        required = FALSE
      ),
      score = ellmer::type_string(
        "Optional numeric combo phenotype column. Auto-detected if omitted.",
        required = FALSE
      ),
      method = ellmer::type_string(
        "One of spearman, pearson, or wilcoxon. Defaults to spearman.",
        required = FALSE
      ),
      match_mode = ellmer::type_string(
        "One of exact or contains. Defaults to exact.",
        required = FALSE
      ),
      symmetric = ellmer::type_boolean(
        "Whether treatment1/treatment2 order can be swapped. Defaults to TRUE.",
        required = FALSE
      ),
      feature_match_mode = ellmer::type_string(
        "One of exact or contains. Defaults to exact.",
        required = FALSE
      ),
      summary_stat = ellmer::type_string(
        "Optional molecular summary.stat override.",
        required = FALSE
      ),
      p_adjust_method = ellmer::type_string(
        "p.adjust method. Defaults to BH.",
        required = FALSE
      )
    )
  ),
  pgx_compare_mono_combo_biomarkers = pgx_tool(
    fun = pgx_compare_mono_combo_biomarkers,
    name = "pgx_compare_mono_combo_biomarkers",
    description = paste(
      "Compare bounded molecular feature associations for a selected combo",
      "phenotype against each monotherapy arm."
    ),
    arguments = list(
      dataset = ellmer::type_string(
        "Dataset alias, dataset name, or explicit local .qs/.rds combo PSet path."
      ),
      treatment1 = ellmer::type_string("First treatment ID/name."),
      treatment2 = ellmer::type_string("Second treatment ID/name."),
      mDataType = ellmer::type_string(
        "Molecular profile type such as rna, mutation, or cnv."
      ),
      features = ellmer::type_array(
        "Feature IDs or symbols to compare.",
        items = ellmer::type_string()
      ),
      combo_table = ellmer::type_string(
        "Combo treatmentResponse table name or auto. Defaults to auto.",
        required = FALSE
      ),
      combo_score = ellmer::type_string(
        "Optional combo phenotype column such as Bliss_score or ZIP_delta.",
        required = FALSE
      ),
      mono_table = ellmer::type_string(
        "Mono treatmentResponse table name. Defaults to mono_profiles.",
        required = FALSE
      ),
      mono_metric = ellmer::type_string(
        "Monotherapy metric column. Defaults to aac_recomputed.",
        required = FALSE
      ),
      method = ellmer::type_string(
        "One of spearman, pearson, or wilcoxon. Defaults to spearman.",
        required = FALSE
      ),
      match_mode = ellmer::type_string(
        "One of exact or contains. Defaults to exact.",
        required = FALSE
      ),
      symmetric = ellmer::type_boolean(
        "Whether treatment1/treatment2 order can be swapped. Defaults to TRUE.",
        required = FALSE
      ),
      feature_match_mode = ellmer::type_string(
        "One of exact or contains. Defaults to exact.",
        required = FALSE
      ),
      summary_stat = ellmer::type_string(
        "Optional molecular summary.stat override.",
        required = FALSE
      ),
      p_adjust_method = ellmer::type_string(
        "p.adjust method. Defaults to BH.",
        required = FALSE
      )
    )
  ),
  pgx_detect_assay_mode = pgx_tool(
    fun = pgx_detect_assay_mode,
    name = "pgx_detect_assay_mode",
    description = "Detect whether a bundled PharmacoGx dataset appears monotherapy, combination, mixed, or perturbation-only.",
    arguments = list(
      dataset = ellmer::type_string(
        "One of GDSCsmall, CCLEsmall, or CMAPsmall. Defaults to GDSCsmall.",
        required = FALSE
      )
    )
  ),
  pgx_find_pset_overlap = pgx_tool(
    fun = pgx_find_pset_overlap,
    name = "pgx_find_pset_overlap",
    description = "Find shared samples, treatments, and sample-treatment pairs across bundled PharmacoGx datasets.",
    arguments = list(
      datasets = ellmer::type_array(
        "Dataset names to compare.",
        items = ellmer::type_string(),
        required = FALSE
      ),
      limit = ellmer::type_integer(
        "Maximum number of shared values to return. Defaults to 50.",
        required = FALSE
      )
    )
  ),
  pgx_compare_pset_response = pgx_tool(
    fun = pgx_compare_pset_response,
    name = "pgx_compare_pset_response",
    description = paste(
      "Compare response values for the same drug and metric across bundled or",
      "local PharmacoSets using exact sample IDs."
    ),
    arguments = list(
      datasets = ellmer::type_array(
        "Dataset names, local manifest aliases, or explicit .qs/.rds paths.",
        items = ellmer::type_string(),
        required = FALSE
      ),
      drug = ellmer::type_string("Treatment or drug name."),
      metric = ellmer::type_string(
        "Sensitivity metric. Defaults to auc_recomputed.",
        required = FALSE
      ),
      samples = ellmer::type_array(
        "Optional sample IDs to include.",
        items = ellmer::type_string(),
        required = FALSE
      ),
      limit = ellmer::type_integer(
        "Maximum number of records to return. Defaults to 50.",
        required = FALSE
      )
    )
  ),
  pgx_summarize_sensitivity = pgx_tool(
    fun = pgx_summarize_sensitivity,
    name = "pgx_summarize_sensitivity",
    description = paste(
      "Summarize a PharmacoSet sensitivity matrix into treatment-sample",
      "records using a selected sensitivity measure."
    ),
    arguments = list(
      dataset = ellmer::type_string(
        paste(
          "Bundled dataset, local manifest alias, or explicit .qs/.rds path.",
          "Defaults to GDSCsmall."
        ),
        required = FALSE
      ),
      sensitivity_measure = ellmer::type_string(
        paste(
          "Sensitivity measure such as auc_recomputed or ic50_recomputed.",
          "Unavailable measures return a catalog instead of raising an error."
        ),
        required = FALSE
      ),
      drugs = ellmer::type_array(
        "Optional vector of treatment names to include.",
        items = ellmer::type_string(),
        required = FALSE
      ),
      cell_lines = ellmer::type_array(
        "Optional vector of sample or cell-line names to include.",
        items = ellmer::type_string(),
        required = FALSE
      ),
      limit = ellmer::type_integer(
        "Maximum number of non-missing records to return. Defaults to 50.",
        required = FALSE
      )
    )
  ),
  pgx_top_responders = pgx_tool(
    fun = pgx_top_responders,
    name = "pgx_top_responders",
    description = paste(
      "Rank samples for one drug using a PharmacoGx sensitivity measure.",
      "Use this for simple responder/non-responder demo questions."
    ),
    arguments = list(
      dataset = ellmer::type_string(
        paste(
          "Bundled dataset, local manifest alias, or explicit .qs/.rds path.",
          "Defaults to GDSCsmall."
        ),
        required = FALSE
      ),
      drug = ellmer::type_string("Treatment name to rank."),
      sensitivity_measure = ellmer::type_string(
        "Sensitivity measure such as auc_recomputed or ic50_recomputed.",
        required = FALSE
      ),
      rank_direction = ellmer::type_string(
        "One of auto, lowest, or highest. Defaults to auto.",
        required = FALSE
      ),
      limit = ellmer::type_integer(
        "Maximum number of ranked samples to return. Defaults to 10.",
        required = FALSE
      )
    )
  ),
  pgx_compute_dose_response_metrics = pgx_tool(
    fun = pgx_compute_dose_response_metrics,
    name = "pgx_compute_dose_response_metrics",
    description = paste(
      "Compute PharmacoGx dose-response metrics from concentration and",
      "viability vectors, including AUC, AAC, IC50, and AC50."
    ),
    arguments = list(
      concentration = ellmer::type_array(
        "Numeric concentration vector.",
        items = ellmer::type_number()
      ),
      viability = ellmer::type_array(
        "Numeric viability vector aligned to concentration.",
        items = ellmer::type_number()
      ),
      viability_as_pct = ellmer::type_boolean(
        "Whether viability is provided as 0-100 percentages. Defaults to TRUE.",
        required = FALSE
      ),
      area_type = ellmer::type_string(
        "AUC/AAC area mode: Actual or Fitted. Defaults to Actual.",
        required = FALSE
      )
    )
  ),
  pgx_compute_synergy_reference = pgx_tool(
    fun = pgx_compute_synergy_reference,
    name = "pgx_compute_synergy_reference",
    description = paste(
      "Compute simple Bliss and HSA null reference viabilities for two",
      "monotherapy response vectors."
    ),
    arguments = list(
      viability_1 = ellmer::type_array(
        "First monotherapy viability vector.",
        items = ellmer::type_number()
      ),
      viability_2 = ellmer::type_array(
        "Second monotherapy viability vector.",
        items = ellmer::type_number()
      ),
      hsa_na_rm = ellmer::type_boolean(
        "Whether computeHSA should ignore missing values. Defaults to FALSE.",
        required = FALSE
      )
    )
  ),
  pgx_list_available_psets = pgx_tool(
    fun = pgx_list_available_psets,
    name = "pgx_list_available_psets",
    description = paste(
      "Fetch the remote PharmacoGx table of downloadable PharmacoSets.",
      "Use this before downloading a full external PSet."
    ),
    arguments = list(
      canonical = ellmer::type_boolean(
        "Whether to list only canonical PSets. Defaults to TRUE.",
        required = FALSE
      ),
      limit = ellmer::type_integer(
        "Maximum number of rows to return. Defaults to 25.",
        required = FALSE
      )
    )
  ),
  pgx_download_pset = pgx_tool(
    fun = pgx_download_pset,
    name = "pgx_download_pset",
    description = paste(
      "Download a named PharmacoSet with PharmacoGx::downloadPSet() and",
      "return JSON-safe metadata. Requires confirm_download = TRUE."
    ),
    arguments = list(
      name = ellmer::type_string(
        "PSet Name exactly as returned by pgx_list_available_psets()."
      ),
      save_dir = ellmer::type_string(
        "Directory where the downloaded PSet should be saved. Defaults to tempdir().",
        required = FALSE
      ),
      pset_file_name = ellmer::type_string(
        "Optional file name for the downloaded .rds file.",
        required = FALSE
      ),
      timeout = ellmer::type_number(
        "Download timeout in seconds. Defaults to 600.",
        required = FALSE
      ),
      confirm_download = ellmer::type_boolean(
        "Must be TRUE to confirm the user approved this external download.",
        required = FALSE
      )
    )
  ),
  pgx_session_info = pgx_tool(
    fun = pgx_session_info,
    name = "pgx_session_info",
    description = "Return R and package version metadata for the PharmacoGx MCP demo.",
    arguments = list()
  )
)

res
