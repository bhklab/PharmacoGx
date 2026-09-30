# PharmacoGx MCP execution

Use this route only when the host actually exposes the typed `pgx_*` tools.
Prefer these tools to terminal or generic-R substitutions and keep one
persistent MCP session for multi-step work.

## Start and discovery

1. Call `pgx_session_info` once to confirm the active package source, tool
   registry, local discovery roots, and whether the server needs a restart.
2. Use `pgx_list_example_datasets` for fixtures or `pgx_list_local_psets` for
   real local PSets.
3. Call `pgx_check_analysis_requirements` once before choosing among candidate
   PSets. Inspect identifier matches, missing requirements, closest checked
   candidates, and datasets that could not be loaded.
4. Use `pgx_list_available_psets` only to list downloadable data. Its output
   does not establish treatment, sample, endpoint, or molecular coverage.
   Call `pgx_download_pset` only after approval and with
   `confirm_download = TRUE`.
5. Call `pgx_list_entities` and `pgx_list_available_covariates` before assuming
   entity names, metrics, or annotation fields.

Do not repeat unchanged discovery calls. If a typed workflow is unavailable,
report the missing capability instead of silently switching to a shell or
generic R.

## Tool routing

### Metadata and filtering

- `pgx_get_sample_metadata`
- `pgx_get_treatment_metadata`
- `pgx_filter_samples`

Use exact covariate names returned by `pgx_list_available_covariates`.

### Drug response and curves

- `pgx_summarize_sensitivity`
- `pgx_top_responders`
- `pgx_compare_metrics`
- `pgx_get_dose_response_points`
- `pgx_plot_dose_response`
- `pgx_compute_dose_response_metrics`

### Biomarker analysis

- `pgx_get_molecular_profile`
- `pgx_suggest_biomarker_candidates` for hypothesis generation only
- `pgx_association_test`
- `pgx_multi_pset_association_test`
- `pgx_find_pset_overlap`
- `pgx_compare_pset_response`

Use `pgx_multi_pset_association_test` for independent per-PSet analyses and
consensus summarization, not raw-data pooling.

### Combination analysis

- `pgx_detect_assay_mode`
- `pgx_list_treatment_response_tables`
- `pgx_rank_synergy_combinations`
- `pgx_get_combo_response_records`
- `pgx_plot_synergy_heatmap`
- `pgx_combo_biomarker_association`
- `pgx_compare_mono_combo_biomarkers`
- `pgx_compute_synergy_reference` when aligned monotherapy viability vectors
  are available

### PSet curation

- `pgx_pset_curation_questions`
- `pgx_validate_pset_inputs`

Do not attempt construction until the inputs pass validation and the user asks
for construction.
