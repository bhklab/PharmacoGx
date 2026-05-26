# PharmacoGx MCP tool definitions.
#
# This file is sourced by mcptools::mcp_server(tools = ...). It must evaluate to
# a named list of ellmer::tool() objects.

pgx_supported_datasets <- c("GDSCsmall", "CCLEsmall", "CMAPsmall")

pgx_load_dataset <- function(dataset) {
  if (!dataset %in% pgx_supported_datasets) {
    stop(
      "Unsupported dataset '",
      dataset,
      "'. Use one of: ",
      paste(pgx_supported_datasets, collapse = ", ")
    )
  }

  env <- new.env(parent = emptyenv())
  utils::data(list = dataset, package = "PharmacoGx", envir = env)
  env[[dataset]]
}

pgx_dataset_catalog <- function() {
  data.frame(
    dataset = c("GDSCsmall", "CCLEsmall", "CMAPsmall", "HDAC_genes"),
    object_type = c("PharmacoSet", "PharmacoSet", "PharmacoSet", "data.frame"),
    analysis_mode = c(
      "drug sensitivity toy dataset",
      "drug sensitivity toy dataset",
      "drug perturbation toy dataset",
      "example HDAC inhibitor gene signature"
    ),
    recommended_demo_use = c(
      "summarize sensitivity profiles and rank responders",
      "inspect cross-study sensitivity and molecular profile examples",
      "connectivity and perturbation signature examples",
      "signature lookup with CMAPsmall"
    ),
    stringsAsFactors = FALSE
  )
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
  x[nzchar(x)]
}

pgx_is_blank <- function(x) {
  is.null(x) || length(x) == 0 || is.na(x[[1]]) || !nzchar(x[[1]])
}

pgx_matrix_to_records <- function(x, value_name = "value", limit = 50) {
  records <- as.data.frame(as.table(x), stringsAsFactors = FALSE)
  names(records) <- c("treatment", "sample", value_name)
  records <- records[!is.na(records[[value_name]]), , drop = FALSE]
  pgx_limit_records(records, limit)
}

pgx_table_with_id <- function(x, id_name) {
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  data.frame(
    stats::setNames(list(rownames(x)), id_name),
    x,
    row.names = NULL,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
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
          tolower(as.character(searchable[[field]])) == tolower(query)
        })
      )
    } else {
      field_match <- Reduce(
        `|`,
        lapply(search_fields, function(field) {
          grepl(
            tolower(query),
            tolower(as.character(searchable[[field]])),
            fixed = TRUE
          )
        })
      )
    }
    matched <- matched | field_match
  }

  which(matched)
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

  values <- as.character(values)
  match_mode <- pgx_resolve_match_mode(match_mode)

  Reduce(
    `|`,
    lapply(expected, function(value) {
      if (identical(match_mode, "exact")) {
        tolower(values) == tolower(value)
      } else {
        grepl(tolower(value), tolower(values), fixed = TRUE)
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

pgx_sensitivity_records <- function(pset, metric, drug = NULL) {
  args <- list(
    object = pset,
    sensitivity.measure = metric,
    verbose = FALSE
  )

  if (!pgx_is_blank(drug)) {
    args$drugs <- drug
  }

  summary <- do.call(PharmacoGx::summarizeSensitivityProfiles, args)
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

pgx_list_example_datasets <- function() {
  list(
    datasets = pgx_dataset_catalog(),
    note = paste(
      "These are bundled demonstration datasets. Use downloadPSet() only when",
      "an explicit large external dataset workflow is intended."
    )
  )
}

pgx_list_entities <- function(
  dataset = "GDSCsmall",
  entity = "treatments",
  limit = 25
) {
  entity <- match.arg(
    entity,
    c("treatments", "samples", "sensitivity_measures", "datasets")
  )
  limit <- max(1L, as.integer(limit))

  if (identical(entity, "datasets")) {
    return(pgx_list_example_datasets())
  }

  pset <- pgx_load_dataset(dataset)
  values <- switch(
    entity,
    "treatments" = PharmacoGx::treatmentNames(pset),
    "samples" = PharmacoGx::sampleNames(pset),
    "sensitivity_measures" = PharmacoGx::sensitivityMeasures(pset)
  )

  list(
    dataset = dataset,
    entity = entity,
    total = length(values),
    returned = min(length(values), limit),
    values = pgx_trim(values, limit)
  )
}

pgx_list_available_covariates <- function(dataset = "GDSCsmall") {
  pset <- pgx_load_dataset(dataset)
  molecular_profiles <- tryCatch(
    PharmacoGx::mDataNames(pset),
    error = function(e) character()
  )

  feature_fields <- lapply(molecular_profiles, function(mDataType) {
    tryCatch(
      colnames(PharmacoGx::featureInfo(pset, mDataType)),
      error = function(e) character()
    )
  })
  names(feature_fields) <- molecular_profiles

  list(
    dataset = dataset,
    sample_fields = colnames(PharmacoGx::sampleInfo(pset)),
    treatment_fields = colnames(PharmacoGx::treatmentInfo(pset)),
    sensitivity_info_fields = colnames(PharmacoGx::sensitivityInfo(pset)),
    sensitivity_measures = PharmacoGx::sensitivityMeasures(pset),
    molecular_profiles = molecular_profiles,
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
  info <- PharmacoGx::sampleInfo(pset)
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
  info <- PharmacoGx::treatmentInfo(pset)
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
      "pubchem"
    )
  )
  matched_rows <- pgx_find_rows(info, treatments, match_mode, search_fields)
  result <- pgx_select_fields(
    info[matched_rows, , drop = FALSE],
    "treatment",
    fields
  )
  limit <- max(1L, as.integer(limit))

  list(
    dataset = dataset,
    query = pgx_clean_vector(treatments),
    match_mode = pgx_resolve_match_mode(match_mode),
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
  info <- PharmacoGx::sampleInfo(pset)
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
  available <- PharmacoGx::sensitivityMeasures(pset)
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
  records_by_metric <- lapply(metrics, function(metric) {
    pgx_sensitivity_records(pset, metric, drug)
  })
  names(records_by_metric) <- metrics

  joined <- Reduce(
    function(x, y) merge(x, y, by = c("treatment", "sample"), all = FALSE),
    records_by_metric
  )

  pairwise <- do.call(
    rbind,
    lapply(utils::combn(metrics, 2, simplify = FALSE), function(pair) {
      pair_values <- joined[, pair, drop = FALSE]
      complete <- stats::complete.cases(pair_values) &
        rowSums(!is.finite(as.matrix(pair_values))) == 0
      data.frame(
        metric_1 = pair[[1]],
        metric_2 = pair[[2]],
        complete_pairs = sum(complete),
        spearman = if (sum(complete) >= 2) {
          stats::cor(
            joined[[pair[[1]]]][complete],
            joined[[pair[[2]]]][complete],
            method = "spearman"
          )
        } else {
          NA_real_
        },
        pearson = if (sum(complete) >= 2) {
          stats::cor(
            joined[[pair[[1]]]][complete],
            joined[[pair[[2]]]][complete],
            method = "pearson"
          )
        } else {
          NA_real_
        },
        stringsAsFactors = FALSE
      )
    })
  )

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
    correlations = pairwise,
    overlap = overlap,
    top_sensitive_records = top_sensitive_records,
    top_resistant_records = top_resistant_records,
    note = paste(
      "This is a descriptive robustness comparison. It is not a formal",
      "statistical validation."
    )
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

  args <- list(
    object = pset,
    sensitivity.measure = sensitivity_measure,
    verbose = FALSE
  )

  if (!is.null(drugs)) {
    args$drugs <- drugs
  }
  if (!is.null(cell_lines)) {
    args$cell.lines <- cell_lines
  }

  summary <- do.call(PharmacoGx::summarizeSensitivityProfiles, args)

  list(
    dataset = dataset,
    sensitivity_measure = sensitivity_measure,
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

  summary <- PharmacoGx::summarizeSensitivityProfiles(
    object = pset,
    sensitivity.measure = sensitivity_measure,
    drugs = drug,
    verbose = FALSE
  )

  if (!drug %in% rownames(summary)) {
    stop(
      "Drug '",
      drug,
      "' was not found in summarized sensitivity results for ",
      dataset,
      "."
    )
  }

  values <- summary[drug, , drop = TRUE]
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
    sensitivity_measure = sensitivity_measure,
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

pgx_compute_synergy_reference <- function(viability_1, viability_2) {
  if (
    !length(viability_1) %in% c(1L, length(viability_2)) &&
      !length(viability_2) %in% c(1L, length(viability_1))
  ) {
    stop("viability_1 and viability_2 must have length 1 or the same length.")
  }

  viability_1 <- as.numeric(viability_1)
  viability_2 <- as.numeric(viability_2)

  data.frame(
    viability_1 = viability_1,
    viability_2 = viability_2,
    bliss_reference = PharmacoGx::computeBliss(viability_1, viability_2),
    hsa_reference = PharmacoGx::computeHSA(viability_1, viability_2)
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
    PharmacoGx::sensitivityMeasures(pset),
    character()
  )
  molecular_profiles <- capture_value(PharmacoGx::mDataNames(pset), character())

  list(
    name = capture_value(PharmacoGx::name(pset), NA_character_),
    class = paste(class(pset), collapse = ", "),
    file_path = file_path,
    dataset_type = capture_value(PharmacoGx::datasetType(pset), NA_character_),
    sample_count = capture_value(
      length(PharmacoGx::sampleNames(pset)),
      NA_integer_
    ),
    treatment_count = capture_value(
      length(PharmacoGx::treatmentNames(pset)),
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
    save_dir <- tempdir()
  }

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

  list(
    downloaded = TRUE,
    requested_pset = name,
    save_dir = normalizePath(save_dir, mustWork = FALSE),
    file_path = normalizePath(file_path, mustWork = FALSE),
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

  list(
    r_version = as.character(getRversion()),
    package_versions = data.frame(
      package = names(versions),
      version = unname(versions),
      installed = !is.na(versions),
      stringsAsFactors = FALSE
    ),
    working_directory = getwd(),
    supported_demo_datasets = pgx_supported_datasets
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
  pgx_list_example_datasets = ellmer::tool(
    fun = pgx_list_example_datasets,
    name = "pgx_list_example_datasets",
    description = paste(
      "List bundled PharmacoGx demo datasets and the type of agentic",
      "analysis each supports."
    ),
    arguments = list()
  ),
  pgx_list_entities = ellmer::tool(
    fun = pgx_list_entities,
    name = "pgx_list_entities",
    description = paste(
      "List treatments, samples, sensitivity measures, or demo datasets",
      "available from a bundled PharmacoGx dataset."
    ),
    arguments = list(
      dataset = ellmer::type_string(
        "One of GDSCsmall, CCLEsmall, or CMAPsmall. Defaults to GDSCsmall.",
        required = FALSE
      ),
      entity = ellmer::type_string(
        paste(
          "Entity to list: treatments, samples, sensitivity_measures,",
          "or datasets. Defaults to treatments."
        ),
        required = FALSE
      ),
      limit = ellmer::type_integer(
        "Maximum number of values to return. Defaults to 25.",
        required = FALSE
      )
    )
  ),
  pgx_list_available_covariates = ellmer::tool(
    fun = pgx_list_available_covariates,
    name = "pgx_list_available_covariates",
    description = paste(
      "List sample, treatment, sensitivity, and molecular covariates",
      "available in a bundled PharmacoGx dataset."
    ),
    arguments = list(
      dataset = ellmer::type_string(
        "One of GDSCsmall, CCLEsmall, or CMAPsmall. Defaults to GDSCsmall.",
        required = FALSE
      )
    )
  ),
  pgx_get_sample_metadata = ellmer::tool(
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
  pgx_get_treatment_metadata = ellmer::tool(
    fun = pgx_get_treatment_metadata,
    name = "pgx_get_treatment_metadata",
    description = paste(
      "Return treatment or drug annotation rows from a bundled PharmacoGx",
      "dataset, including target and mechanism fields when available."
    ),
    arguments = list(
      dataset = ellmer::type_string(
        "One of GDSCsmall, CCLEsmall, or CMAPsmall. Defaults to GDSCsmall.",
        required = FALSE
      ),
      treatments = ellmer::type_array(
        "Treatment or drug IDs/names to look up.",
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
  pgx_filter_samples = ellmer::tool(
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
      filters = ellmer::type_object(
        .description = "Named sample metadata filters, e.g. {\"tissueid\": \"lung\"}.",
        .additional_properties = TRUE
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
  pgx_compare_metrics = ellmer::tool(
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
  pgx_summarize_sensitivity = ellmer::tool(
    fun = pgx_summarize_sensitivity,
    name = "pgx_summarize_sensitivity",
    description = paste(
      "Summarize a PharmacoSet sensitivity matrix into treatment-sample",
      "records using a selected sensitivity measure."
    ),
    arguments = list(
      dataset = ellmer::type_string(
        "One of GDSCsmall or CCLEsmall. Defaults to GDSCsmall.",
        required = FALSE
      ),
      sensitivity_measure = ellmer::type_string(
        "Sensitivity measure such as auc_recomputed, auc_published, or ic50_recomputed.",
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
  pgx_top_responders = ellmer::tool(
    fun = pgx_top_responders,
    name = "pgx_top_responders",
    description = paste(
      "Rank samples for one drug using a PharmacoGx sensitivity measure.",
      "Use this for simple responder/non-responder demo questions."
    ),
    arguments = list(
      dataset = ellmer::type_string(
        "One of GDSCsmall or CCLEsmall. Defaults to GDSCsmall.",
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
  pgx_compute_dose_response_metrics = ellmer::tool(
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
  pgx_compute_synergy_reference = ellmer::tool(
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
      )
    )
  ),
  pgx_list_available_psets = ellmer::tool(
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
  pgx_download_pset = ellmer::tool(
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
  pgx_session_info = ellmer::tool(
    fun = pgx_session_info,
    name = "pgx_session_info",
    description = "Return R and package version metadata for the PharmacoGx MCP demo.",
    arguments = list()
  )
)

res
