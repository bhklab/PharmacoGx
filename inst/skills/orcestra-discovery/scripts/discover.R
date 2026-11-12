# Read-only discovery adapter. Source this file, then call pgx_discover().
pgx_discovery_request <- function(kind, value = NULL) {
  kind <- match.arg(
    kind,
    c("datasets", "dataset_counts", "compound", "gene", "psets")
  )
  if (kind %in% c("dataset_counts", "compound", "gene")) {
    if (
      !is.character(value) ||
        length(value) != 1L ||
        is.na(value) ||
        !nzchar(trimws(value)) ||
        nchar(value) > 200L
    ) {
      stop("Supply one nonempty entity name of at most 200 characters.")
    }
  }
  queries <- c(
    datasets = "query { datasets { id name } }",
    dataset_counts = paste(
      "query Lookup($name: String) { dataset(datasetName: $name) {",
      "id name cell_count compound_tested_count tissue_tested_count experiment_count } }"
    ),
    compound = paste(
      "query Lookup($name: String) { compound(compoundName: $name) {",
      "compound { id uid name } } }"
    ),
    gene = "query Lookup($name: String) { gene(geneName: $name) { id name annotation { symbol } } }"
  )
  list(
    operation = kind,
    backend = if (kind == "psets") {
      "orcestra-pset-catalog"
    } else {
      "pharmacodb-graphql"
    },
    query = if (kind == "psets") NULL else unname(queries[[kind]]),
    variables = if (kind %in% c("datasets", "psets")) {
      list()
    } else {
      list(name = value)
    }
  )
}

pgx_discovery_payload_status <- function(data, operation = NULL) {
  if (!is.list(data) || is.null(names(data))) {
    return("invalid_shape")
  }
  fields <- c(
    datasets = "datasets",
    dataset_counts = "dataset",
    compound = "compound",
    gene = "gene"
  )
  if (is.null(operation)) {
    matches <- names(fields)[fields %in% names(data)]
    if (length(matches) != 1L) {
      return("missing_field")
    }
    operation <- matches[[1L]]
  }
  if (!operation %in% names(fields)) {
    return("invalid_shape")
  }
  field <- unname(fields[[operation]])
  if (!field %in% names(data)) {
    return("missing_field")
  }
  value <- data[[field]]
  if (is.null(value)) {
    return("empty")
  }
  if (operation == "compound") {
    if (!is.list(value) || is.null(names(value))) {
      return("invalid_shape")
    }
    if (!"compound" %in% names(value)) {
      return("missing_field")
    }
    value <- value$compound
    if (is.null(value)) return("empty")
  }
  collection <- operation %in% c("datasets", "dataset_counts")
  if (collection) {
    if (!is.list(value) || !is.null(names(value))) {
      return("invalid_shape")
    }
    if (!length(value)) return("empty")
  } else {
    value <- list(value)
  }
  required <- switch(
    operation,
    datasets = c("id", "name"),
    dataset_counts = c(
      "id",
      "name",
      "cell_count",
      "compound_tested_count",
      "tissue_tested_count",
      "experiment_count"
    ),
    compound = c("id", "uid", "name"),
    gene = c("id", "name", "annotation")
  )
  valid <- vapply(
    value,
    function(row) {
      is.list(row) &&
        all(required %in% names(row)) &&
        is.character(row$name) &&
        length(row$name) == 1L &&
        !is.na(row$name) &&
        nzchar(row$name) &&
        length(row$id) == 1L &&
        !is.na(row$id)
    },
    logical(1)
  )
  if (all(valid)) "ok" else "invalid_shape"
}

pgx_decode_discovery <- function(
  body,
  http_status,
  graphql = TRUE,
  operation = NULL
) {
  result <- list(
    status = "ok",
    errors = list(),
    data = NULL,
    raw_response = body
  )
  if (http_status < 200L || http_status >= 300L) {
    result$status <- "http_error"
    result$errors <- list(paste("HTTP", http_status))
    return(result)
  }
  parsed <- tryCatch(
    jsonlite::fromJSON(body, simplifyVector = FALSE),
    error = identity
  )
  if (inherits(parsed, "error")) {
    result$status <- "invalid_json"
    result$errors <- list(conditionMessage(parsed))
    return(result)
  }
  if (graphql) {
    if (!is.list(parsed) || is.null(names(parsed))) {
      result$status <- "invalid_shape"
    } else if (length(parsed$errors) > 0L) {
      result$status <- "graphql_error"
      result$errors <- parsed$errors
    } else if (is.null(parsed$data)) {
      result$status <- "missing_data"
    } else {
      result$status <- pgx_discovery_payload_status(parsed$data, operation)
      if (result$status == "ok") result$data <- parsed$data
    }
  } else if (
    !is.list(parsed) ||
      !is.null(names(parsed)) ||
      any(
        !vapply(
          parsed,
          function(x) {
            is.list(x) &&
              is.character(x$name) &&
              length(x$name) == 1L
          },
          logical(1)
        )
      )
  ) {
    result$status <- "invalid_shape"
  } else if (!length(parsed)) {
    result$status <- "empty"
  } else {
    result$data <- parsed
  }
  result
}

pgx_discover <- function(kind, value = NULL, endpoint = NULL) {
  if (
    !requireNamespace("curl", quietly = TRUE) ||
      !requireNamespace("jsonlite", quietly = TRUE)
  ) {
    stop("Discovery requires curl and jsonlite in the selected R environment.")
  }
  request <- pgx_discovery_request(kind, value)
  graphql <- request$backend == "pharmacodb-graphql"
  if (is.null(endpoint)) {
    endpoint <- if (graphql) {
      "https://pharmacodb.ca/graphql"
    } else {
      "https://orcestra.ca/api/pset/available"
    }
  }
  if (
    length(endpoint) != 1L ||
      is.na(endpoint) ||
      !grepl("^https://[^/@]+(/|$)", endpoint)
  ) {
    stop("Use one verified HTTPS endpoint without embedded credentials.")
  }
  handle <- curl::new_handle(
    timeout = 30,
    connecttimeout = 10,
    maxfilesize = 2 * 1024^2,
    followlocation = FALSE
  )
  if (graphql) {
    payload <- list(query = request$query)
    if (length(request$variables)) {
      payload$variables <- request$variables
    }
    curl::handle_setopt(
      handle,
      postfields = jsonlite::toJSON(payload, auto_unbox = TRUE)
    )
    curl::handle_setheaders(handle, "Content-Type" = "application/json")
  }
  envelope <- c(
    request,
    list(
      adapter_schema_version = 1L,
      guide_revision = if (graphql) {
        "193e18c4a93d46e1906bb1195ec0d3bac4758cb2"
      } else {
        NULL
      },
      backend_release = NULL,
      endpoint = endpoint,
      retrieved_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
    )
  )
  response <- tryCatch(
    curl::curl_fetch_memory(endpoint, handle),
    error = identity
  )
  if (inherits(response, "error")) {
    return(c(
      envelope,
      list(
        http_status = NULL,
        status = "transport_error",
        errors = list(conditionMessage(response)),
        data = NULL,
        raw_response = NULL
      )
    ))
  }
  if (length(response$content) > 2 * 1024^2) {
    stop("Response exceeds the 2 MiB budget.")
  }
  c(
    envelope,
    list(http_status = response$status_code),
    pgx_decode_discovery(
      rawToChar(response$content),
      response$status_code,
      graphql,
      operation = request$operation
    )
  )
}
