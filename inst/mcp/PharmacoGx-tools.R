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

pgx_safe_filename <- function(...) {
  x <- paste(..., sep = "_")
  x <- gsub("[^A-Za-z0-9_.-]+", "_", x)
  x <- gsub("_+", "_", x)
  x
}

pgx_resolve_experiments <- function(pset, drug, sample) {
  info <- PharmacoGx::sensitivityInfo(pset)
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

  raw <- PharmacoGx::sensitivityRaw(pset)
  if (is.null(dim(raw)) || length(dim(raw)) != 3 || any(dim(raw) == 0)) {
    stop(
      "No sensitivityRaw() dose-response array is available for this dataset."
    )
  }

  info <- PharmacoGx::sensitivityInfo(pset)
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

  annotation <- tryCatch(
    S4Vectors::metadata(PharmacoGx::molecularProfilesSlot(pset)[[
      mDataType
    ]])$annotation,
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
  feature_info <- as.data.frame(
    PharmacoGx::featureInfo(pset, mDataType)
  )
  match_fields <- colnames(feature_info)
  matched <- pgx_find_rows(feature_info, features, match_mode, match_fields)
  if (length(matched) == 0) {
    return(list(
      feature_ids = character(),
      matched_features = data.frame()
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
    matched_features = matched_features
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
  samples <- pgx_clean_vector(samples)
  if (is.null(samples)) {
    samples <- PharmacoGx::sampleNames(pset)
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

  mat <- PharmacoGx::summarizeMolecularProfiles(
    pset,
    mDataType = mDataType,
    features = resolved$feature_ids,
    cell.lines = samples,
    summary.stat = summary_stat,
    verbose = FALSE
  )
  mat <- pgx_plain_matrix(mat)

  list(
    matrix = mat,
    matched_features = resolved$matched_features,
    summary_stat = summary_stat
  )
}

pgx_matrix_records <- function(mat, row_name = "feature", col_name = "sample") {
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
  info <- PharmacoGx::sensitivityInfo(pset)
  if (!all(c("sampleid", "treatmentid") %in% colnames(info))) {
    return(data.frame(sample = character(), treatment = character()))
  }

  unique(data.frame(
    sample = info$sampleid,
    treatment = info$treatmentid,
    stringsAsFactors = FALSE
  ))
}

pgx_load_many_datasets <- function(datasets) {
  datasets <- pgx_clean_vector(datasets)
  if (is.null(datasets) || length(datasets) < 2) {
    stop("At least two datasets are required.")
  }
  stats::setNames(lapply(datasets, pgx_load_dataset), datasets)
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
  limit_samples = 25
) {
  pset <- pgx_load_dataset(dataset)
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
    summary_stat = molecular$summary_stat,
    matched_features = molecular$matched_features,
    dimensions = stats::setNames(
      as.integer(dim(molecular$matrix)),
      c("features", "samples")
    ),
    records = pgx_matrix_records(molecular$matrix)
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
  response <- PharmacoGx::summarizeSensitivityProfiles(
    pset,
    sensitivity.measure = metric,
    drugs = drug,
    verbose = FALSE
  )
  if (!drug %in% rownames(response)) {
    stop("Drug '", drug, "' was not found in sensitivity summaries.")
  }

  response_values <- response[drug, , drop = TRUE]
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

  results <- do.call(
    rbind,
    lapply(rownames(mat), function(feature_id) {
      common_samples <- intersect(names(response_values), colnames(mat))
      x <- as.vector(mat[feature_id, common_samples, drop = TRUE])
      y <- as.numeric(response_values[common_samples])
      complete <- stats::complete.cases(x, y) & is.finite(y)
      x <- x[complete]
      y <- y[complete]

      if (length(y) < 3 || length(unique(x)) < 2) {
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
        test <- suppressWarnings(stats::cor.test(x, y, method = method))
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
        effect <- diff(tapply(y, group, median, na.rm = TRUE))
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

  list(
    dataset = dataset,
    drug = drug,
    metric = metric,
    mDataType = mDataType,
    matched_features = molecular$matched_features,
    results = results,
    note = paste(
      "Exploratory association only. Bundled small datasets are demo fixtures",
      "and are not biologically powered for biomarker claims."
    )
  )
}

pgx_detect_assay_mode <- function(dataset = "GDSCsmall") {
  pset <- pgx_load_dataset(dataset)
  info <- PharmacoGx::sensitivityInfo(pset)
  treatment_names <- tryCatch(
    PharmacoGx::treatmentNames(pset),
    error = function(e) character()
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
  has_combo <- has_combo_cols && combo_rows > 0 || combo_names > 0
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
    note = if (has_combo) {
      "Combination-like fields were detected; inspect raw data before synergy analysis."
    } else {
      "No combination-dose fields were detected in this bundled dataset."
    }
  )
}

pgx_find_pset_overlap <- function(
  datasets = c("GDSCsmall", "CCLEsmall"),
  limit = 50
) {
  psets <- pgx_load_many_datasets(datasets)
  limit <- max(1L, as.integer(limit))

  sample_sets <- lapply(psets, PharmacoGx::sampleNames)
  treatment_sets <- lapply(psets, PharmacoGx::treatmentNames)
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
      if (!metric %in% PharmacoGx::sensitivityMeasures(psets[[dataset]])) {
        return(data.frame())
      }
      mat <- PharmacoGx::summarizeSensitivityProfiles(
        psets[[dataset]],
        sensitivity.measure = metric,
        drugs = drug,
        verbose = FALSE
      )
      if (!drug %in% rownames(mat)) {
        return(data.frame())
      }
      values <- mat[drug, , drop = TRUE]
      if (!is.null(samples)) {
        values <- values[names(values) %in% samples]
      }
      data.frame(
        dataset = dataset,
        sample = names(values),
        treatment = drug,
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

  data.frame(
    viability_1 = viability_1,
    viability_2 = viability_2,
    bliss_reference = PharmacoGx::computeBliss(viability_1, viability_2),
    hsa_reference = PharmacoGx::computeHSA(
      viability_1,
      viability_2,
      na.rm = isTRUE(hsa_na_rm)
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
  pgx_get_dose_response_points = ellmer::tool(
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
  pgx_plot_dose_response = ellmer::tool(
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
  pgx_pset_curation_questions = ellmer::tool(
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
  pgx_validate_pset_inputs = ellmer::tool(
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
  pgx_get_molecular_profile = ellmer::tool(
    fun = pgx_get_molecular_profile,
    name = "pgx_get_molecular_profile",
    description = paste(
      "Return summarized molecular profile values for selected features and",
      "samples from a bundled PharmacoGx dataset."
    ),
    arguments = list(
      dataset = ellmer::type_string(
        "One of GDSCsmall or CCLEsmall. Defaults to GDSCsmall.",
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
      )
    )
  ),
  pgx_association_test = ellmer::tool(
    fun = pgx_association_test,
    name = "pgx_association_test",
    description = paste(
      "Run exploratory feature-response associations for one drug using",
      "summarized molecular profiles and sensitivity metrics."
    ),
    arguments = list(
      dataset = ellmer::type_string(
        "One of GDSCsmall or CCLEsmall. Defaults to GDSCsmall.",
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
  pgx_detect_assay_mode = ellmer::tool(
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
  pgx_find_pset_overlap = ellmer::tool(
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
  pgx_compare_pset_response = ellmer::tool(
    fun = pgx_compare_pset_response,
    name = "pgx_compare_pset_response",
    description = paste(
      "Compare response values for the same drug and metric across bundled",
      "PharmacoGx datasets using exact sample IDs."
    ),
    arguments = list(
      datasets = ellmer::type_array(
        "Dataset names to compare.",
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
      ),
      hsa_na_rm = ellmer::type_boolean(
        "Whether computeHSA should ignore missing values. Defaults to FALSE.",
        required = FALSE
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
