test_that("MCP registry exposes the curated PharmacoGx surface", {
  skip_if_not_installed("ellmer")
  tool_file <- system.file(
    "mcp",
    "PharmacoGx-tools.R",
    package = "PharmacoGx"
  )
  expect_true(file.exists(tool_file))

  tool_environment <- new.env(parent = globalenv())
  tools <- source(tool_file, local = tool_environment)$value

  expect_length(tools, 33)
  expect_setequal(
    names(tools),
    c(
      "pgx_list_example_datasets",
      "pgx_list_local_psets",
      "pgx_check_analysis_requirements",
      "pgx_list_entities",
      "pgx_list_available_covariates",
      "pgx_get_sample_metadata",
      "pgx_get_treatment_metadata",
      "pgx_filter_samples",
      "pgx_compare_metrics",
      "pgx_get_dose_response_points",
      "pgx_plot_dose_response",
      "pgx_pset_curation_questions",
      "pgx_validate_pset_inputs",
      "pgx_get_molecular_profile",
      "pgx_association_test",
      "pgx_suggest_biomarker_candidates",
      "pgx_multi_pset_association_test",
      "pgx_list_treatment_response_tables",
      "pgx_rank_synergy_combinations",
      "pgx_get_combo_response_records",
      "pgx_plot_synergy_heatmap",
      "pgx_combo_biomarker_association",
      "pgx_compare_mono_combo_biomarkers",
      "pgx_detect_assay_mode",
      "pgx_find_pset_overlap",
      "pgx_compare_pset_response",
      "pgx_summarize_sensitivity",
      "pgx_top_responders",
      "pgx_compute_dose_response_metrics",
      "pgx_compute_synergy_reference",
      "pgx_list_available_psets",
      "pgx_download_pset",
      "pgx_session_info"
    )
  )
})

test_that("MCP wrappers return one JSON content block and structured errors", {
  skip_if_not_installed("ellmer")
  tool_file <- system.file(
    "mcp",
    "PharmacoGx-tools.R",
    package = "PharmacoGx"
  )
  tool_environment <- new.env(parent = globalenv())
  source(tool_file, local = tool_environment)

  result <- tool_environment$pgx_as_mcp_result(list(ok = TRUE, value = 1))
  expect_s7_class(result, ellmer::ContentToolResult)
  expect_length(result@value, 1)
  expect_s7_class(result@value[[1]], ellmer::ContentText)
  expect_equal(
    jsonlite::fromJSON(result@value[[1]]@text),
    list(ok = TRUE, value = 1L)
  )

  error_result <- tool_environment$pgx_execute_logged(
    function() stop("Unsupported dataset 'missing'"),
    "pgx_test",
    list(),
    environment()
  )
  expect_false(error_result$ok)
  expect_equal(error_result$error$code, "DATASET_NOT_FOUND")
})

test_that("local PSet roots are discovered without a manifest", {
  skip_if_not_installed("ellmer")
  tool_file <- system.file(
    "mcp",
    "PharmacoGx-tools.R",
    package = "PharmacoGx"
  )
  tool_environment <- new.env(parent = globalenv())
  source(tool_file, local = tool_environment)

  root <- tempfile("pgx-local-psets-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  saveRDS(list(example = TRUE), file.path(root, "PSet_PRISM.rds"))
  saveRDS(list(not_a_pset = TRUE), file.path(root, "normalized_data.rds"))
  catalog <- tool_environment$pgx_local_pset_catalog(
    manifest_path = file.path(root, "missing.csv"),
    roots = root,
    discover = TRUE,
    include_missing = FALSE
  )

  expect_equal(catalog$name, "prism")
  expect_equal(catalog$source, "discovered")
  expect_true(catalog$exists)
})

test_that("entity search supports query and pagination", {
  skip_if_not_installed("ellmer")
  tool_file <- system.file(
    "mcp",
    "PharmacoGx-tools.R",
    package = "PharmacoGx"
  )
  tool_environment <- new.env(parent = globalenv())
  source(tool_file, local = tool_environment)

  result <- tool_environment$pgx_list_entities(
    "GDSCsmall",
    entity = "treatments",
    query = "gemcitabine",
    match_mode = "exact",
    offset = 0,
    limit = 1
  )
  expect_equal(result$total, 1)
  expect_equal(result$values, "Gemcitabine")
  expect_null(result$next_offset)
})

test_that("invalid treatment text does not break metadata lookup", {
  skip_if_not_installed("ellmer")
  tool_file <- system.file(
    "mcp",
    "PharmacoGx-tools.R",
    package = "PharmacoGx"
  )
  tool_environment <- new.env(parent = globalenv())
  source(tool_file, local = tool_environment)

  result <- tool_environment$pgx_get_treatment_metadata(
    "CCLEsmall",
    treatments = "Camptosar",
    match_mode = "contains",
    limit = 5
  )
  expect_equal(result$total_matches, 1)
  expect_match(result$metadata$matched_fields, "brand", ignore.case = TRUE)
})

test_that("analysis preflight reports unmet combined requirements", {
  skip_if_not_installed("ellmer")
  tool_file <- system.file(
    "mcp",
    "PharmacoGx-tools.R",
    package = "PharmacoGx"
  )
  tool_environment <- new.env(parent = globalenv())
  source(tool_file, local = tool_environment)

  result <- tool_environment$pgx_check_analysis_requirements(
    datasets = c("GDSCsmall", "CCLEsmall"),
    drug = "Gemcitabine",
    sample_query = "pancreas",
    metric = "aac_recomputed",
    mDataType = "rna",
    features = "SLC29A1"
  )
  expect_length(result$compatible_datasets, 0)
  expect_false(any(result$requirements$requirements_met))
  expect_true(nrow(result$best_local_candidates) > 0)
  expect_false(result$remote_catalog_assesses_compatibility)
  expect_match(result$next_step, "identifier matches")
  expect_false(grepl("Use pgx_list_available_psets", result$next_step))
})

test_that("local requirements checks isolate full PSet loading", {
  skip_if_not_installed("ellmer")
  tool_file <- system.file(
    "mcp",
    "PharmacoGx-tools.R",
    package = "PharmacoGx"
  )
  probe_script <- system.file(
    "mcp",
    "probe-local-pset.R",
    package = "PharmacoGx"
  )
  expect_true(file.exists(probe_script))

  tool_environment <- new.env(parent = globalenv())
  source(tool_file, local = tool_environment)
  data_environment <- new.env(parent = emptyenv())
  utils::data(
    list = "GDSCsmall",
    package = "PharmacoGx",
    envir = data_environment
  )
  local_pset <- tempfile("GDSCsmall-", fileext = ".rds")
  capability_cache <- tempfile("pgx-capability-cache-")
  saveRDS(data_environment$GDSCsmall, local_pset)
  dir.create(capability_cache)
  on.exit(unlink(c(local_pset, capability_cache), recursive = TRUE), add = TRUE)

  env_names <- c(
    "PHARMACOGX_MCP_CAPABILITY_CACHE_DIR",
    "PHARMACOGX_MCP_PACKAGE_SOURCE",
    "PHARMACOGX_MCP_TOOL_FILE",
    "PHARMACOGX_MCP_DIR"
  )
  old_env <- Sys.getenv(env_names, unset = NA_character_)
  on.exit(
    {
      for (name in env_names) {
        if (is.na(old_env[[name]])) {
          Sys.unsetenv(name)
        } else {
          do.call(Sys.setenv, stats::setNames(list(old_env[[name]]), name))
        }
      }
    },
    add = TRUE
  )
  Sys.setenv(
    PHARMACOGX_MCP_CAPABILITY_CACHE_DIR = capability_cache,
    PHARMACOGX_MCP_PACKAGE_SOURCE = "installed",
    PHARMACOGX_MCP_TOOL_FILE = tool_file,
    PHARMACOGX_MCP_DIR = dirname(tool_file)
  )

  first <- tool_environment$pgx_check_analysis_requirements(
    datasets = local_pset,
    drug = "Doxorubicin",
    include_local = FALSE,
    include_bundled = FALSE
  )
  expect_true(first$requirements$loadable)
  expect_match(first$requirements$check_mode, "^isolated_worker_generated$")
  expect_false(
    paste0("local:", normalizePath(local_pset)) %in%
      names(tool_environment$pgx_dataset_cache)
  )

  second <- tool_environment$pgx_check_analysis_requirements(
    datasets = local_pset,
    drug = "Doxorubicin",
    include_local = FALSE,
    include_bundled = FALSE
  )
  expect_match(second$requirements$check_mode, "^isolated_worker_memory$")
})

test_that("identifier and metric aliases resolve without losing stored IDs", {
  skip_if_not_installed("ellmer")
  tool_file <- system.file(
    "mcp",
    "PharmacoGx-tools.R",
    package = "PharmacoGx"
  )
  tool_environment <- new.env(parent = globalenv())
  source(tool_file, local = tool_environment)

  annotations <- data.frame(
    treatmentid = c("Gemcitabine", "Example drug"),
    cid = c(356653, 1234),
    gene_id = c("ENSG00000112759", "ENSG00000123456.2"),
    aliases = c("dFdC /// Gemzar", "alias-a;alias-b"),
    row.names = c("drug-1", "drug-2")
  )
  pubchem <- tool_environment$pgx_resolve_annotation_queries(
    annotations,
    "PubChem CID: 356653"
  )
  ensembl <- tool_environment$pgx_resolve_annotation_queries(
    annotations,
    "ENSG00000123456.9"
  )
  delimited <- tool_environment$pgx_resolve_annotation_queries(
    annotations,
    "Gemzar"
  )
  metric <- tool_environment$pgx_resolve_metric(
    "auc_recomputed",
    c("aac_recomputed", "ic50_recomputed")
  )

  expect_equal(pubchem$rows, 1)
  expect_equal(ensembl$rows, 2)
  expect_equal(delimited$rows, 1)
  expect_equal(metric$resolved, "aac_recomputed")
  expect_equal(metric$match_type, "legacy_alias")
  expect_true(
    "pancreas" %in% tool_environment$pgx_sample_context_variants("PDAC")
  )
})

test_that("legacy local combo PSets can be classified when available", {
  skip_if_not_installed("ellmer")
  tool_file <- system.file(
    "mcp",
    "PharmacoGx-tools.R",
    package = "PharmacoGx"
  )
  tool_environment <- new.env(parent = globalenv())
  source(tool_file, local = tool_environment)

  local_pset <- test_path(
    "..",
    "..",
    "local_data",
    "NCI_ALMANAC_2017.rds"
  )
  skip_if_not(file.exists(local_pset))
  result <- tool_environment$pgx_detect_assay_mode(local_pset)
  expect_true(result$assay_mode %in% c("combination", "mixed"))
  expect_gt(result$combo_rows, 0)
})

test_that("downloads remain confirmation gated", {
  skip_if_not_installed("ellmer")
  tool_file <- system.file(
    "mcp",
    "PharmacoGx-tools.R",
    package = "PharmacoGx"
  )
  tool_environment <- new.env(parent = globalenv())
  source(tool_file, local = tool_environment)

  result <- tool_environment$pgx_download_pset(
    "example",
    confirm_download = FALSE
  )
  expect_false(result$downloaded)
  expect_true(result$confirmation_required)
})

test_that("downloads report discovery accurately for user-selected directories", {
  skip_if_not_installed("ellmer")
  tool_environment <- new.env(parent = globalenv())
  source(
    system.file("mcp", "PharmacoGx-tools.R", package = "PharmacoGx"),
    local = tool_environment
  )
  destination <- tempfile("download-test-")
  dir.create(destination)
  withr::defer(unlink(destination, recursive = TRUE))
  local_mocked_bindings(
    downloadPSet = function(name, saveDir, pSetFileName, ...) {
      saveRDS(list(test = TRUE), file.path(saveDir, pSetFileName))
      list(test = TRUE)
    },
    .package = "PharmacoGx"
  )
  tool_environment$pgx_pset_summary <- function(...) list(test = TRUE)
  tool_environment$pgx_default_local_pset_roots <- function() character()
  result <- tool_environment$pgx_download_pset(
    "test",
    save_dir = destination,
    pset_file_name = "PSet_test.rds",
    confirm_download = TRUE
  )
  expect_true(result$downloaded)
  expect_true(file.exists(result$file_path))
  expect_false(result$automatically_discoverable)
  tool_environment$pgx_default_local_pset_roots <- function() destination
  result <- tool_environment$pgx_download_pset(
    "test",
    save_dir = destination,
    pset_file_name = "PSet_test.rds",
    confirm_download = TRUE
  )
  expect_true(result$automatically_discoverable)
})

test_that("host-neutral development cases validate", {
  skip_if_not_installed("yaml")
  validator <- system.file(
    "agentic-eval",
    "scripts",
    "validate-cases.R",
    package = "PharmacoGx"
  )
  manifest <- system.file(
    "agentic-eval",
    "development-cases.yaml",
    package = "PharmacoGx"
  )
  expect_true(file.exists(validator))
  expect_true(file.exists(manifest))
  output <- system2(
    file.path(R.home("bin"), "Rscript"),
    c("--vanilla", shQuote(validator), shQuote(manifest)),
    stdout = TRUE,
    stderr = TRUE
  )
  expect_equal(
    attr(output, "status"),
    NULL,
    info = paste(output, collapse = "\n")
  )
  result <- jsonlite::fromJSON(paste(output, collapse = "\n"))
  expect_true(result$valid)
  expect_equal(result$case_count, 10)
})
