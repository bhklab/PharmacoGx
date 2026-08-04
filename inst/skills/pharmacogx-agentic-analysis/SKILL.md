---
name: pharmacogx-agentic-analysis
description: Use PharmacoGx for preclinical pharmacogenomic analysis with MCP-backed R workflows.
---

## Purpose

Use this skill when a user asks an agent to perform, explain, or plan
preclinical pharmacogenomic analysis with PharmacoGx.

Prefer typed PharmacoGx MCP tools over arbitrary R execution. Current tools:

- `pgx_list_example_datasets`
- `pgx_list_local_psets`
- `pgx_list_entities`
- `pgx_list_available_covariates`
- `pgx_get_sample_metadata`
- `pgx_get_treatment_metadata`
- `pgx_filter_samples`
- `pgx_compare_metrics`
- `pgx_get_dose_response_points`
- `pgx_plot_dose_response`
- `pgx_summarize_sensitivity`
- `pgx_top_responders`
- `pgx_compute_dose_response_metrics`
- `pgx_compute_synergy_reference`
- `pgx_pset_curation_questions`
- `pgx_validate_pset_inputs`
- `pgx_get_molecular_profile`
- `pgx_association_test`
- `pgx_suggest_biomarker_candidates`
- `pgx_multi_pset_association_test`
- `pgx_detect_assay_mode`
- `pgx_list_treatment_response_tables`
- `pgx_rank_synergy_combinations`
- `pgx_get_combo_response_records`
- `pgx_plot_synergy_heatmap`
- `pgx_combo_biomarker_association`
- `pgx_compare_mono_combo_biomarkers`
- `pgx_find_pset_overlap`
- `pgx_compare_pset_response`
- `pgx_list_available_psets`
- `pgx_download_pset`
- `pgx_session_info`

## General Workflow

1. Start by identifying the dataset, PharmacoSet, sample, treatment, metric,
   molecular profile, and metadata fields needed for the question.
2. Use bundled toy datasets (`GDSCsmall`, `CCLEsmall`, `CMAPsmall`) only for
   examples, testing, and lightweight demonstrations. For real analysis, accept
   local manifest aliases such as `nci_almanac`, `gdsc2_matrix`, and
   `gdsc2_anchor`, or explicit local `.qs` or `.rds` PharmacoSet file paths.
   Clearly say whether results come from toy data or local PSet-backed data.
   Keep large local PSets, especially `nci_almanac`, in one long-lived session
   across multi-step workflows instead of repeatedly cold-loading them.
3. Use `pgx_list_available_psets` for larger downloadable PharmacoSets only
   when the user asks for external data. Do not call `pgx_download_pset` unless
   the user explicitly approves runtime, storage, and network use.
4. Call `pgx_list_available_covariates` before metadata lookup, filtering, or
   molecular-profile analysis because annotation fields differ by PharmacoSet.
5. Keep outputs bounded. Prefer ranked, filtered, summarized, or plotted
   results over full matrices.
6. Include `pgx_session_info` when summarizing a reproducible run.

## Dataset Inspection

Do not assume all PharmacoSets contain the same samples, treatments, metrics,
annotations, molecular data, fitted curves, or raw dose-response data.

Use:

- `pgx_list_example_datasets` when no dataset is specified.
- `pgx_list_local_psets` before real local PSet-backed demos.
- `pgx_list_entities` before assuming a sample, drug, or metric exists.
- `pgx_list_available_covariates` before using sample, treatment, or molecular
  fields.
- `pgx_detect_assay_mode` before any combination or synergy analysis.
- `pgx_list_treatment_response_tables` when a real combo PSet may store useful
  endpoints in treatmentResponse `profiles`, `synergy`, or `raw` tables.

If the requested dataset, sample, drug, metric, molecular profile, or metadata
field is unavailable, state what is missing and ask for a valid option.

## Metadata Lookup And Filtering

Use metadata tools when users ask about cell lines, tissues, subtypes,
lineages, drugs, targets, mechanisms, or sample/treatment annotations.

Use:

- `pgx_get_sample_metadata` for sample or cell-line annotations.
- `pgx_get_treatment_metadata` for drug or treatment annotations.
- `pgx_filter_samples` to subset samples by tissue, subtype, origin, source, or
  other listed metadata fields.

Use exact covariate names returned by `pgx_list_available_covariates`. Do not
assume fields such as `disease`, `histology`, `lineage`, `age`, or `sex` exist.

Do not call tissue, lineage, or site-level metadata patterns molecular
biomarkers. If only sample metadata were used, describe the result as a
metadata association.

## Drug Sensitivity And Responder Analysis

Use:

- `pgx_top_responders` for most sensitive or resistant samples for a drug.
- `pgx_summarize_sensitivity` for bounded treatment-sample metric tables.
- `pgx_compare_metrics` before claiming a responder pattern is robust.

Explain metric direction:

- Lower AUC usually indicates stronger sensitivity.
- Lower IC50 usually indicates stronger sensitivity when IC50 is reached or
  estimated reliably.
- Higher AAC-like or activity-area metrics usually indicate stronger
  sensitivity.

If metric correlations have unexpected direction, or top-N overlap is
meaningless because N is too large, avoid claiming robustness.

For small sample sizes or singleton metadata groups, describe results as
exploratory.

## Dose-Response Analysis

When users ask about dose-response curves, fitted response metrics, AUC, AAC,
IC50, AC50, Einf, or curve quality, first determine whether they provided raw
vectors or want PSet-backed data.

Use:

- `pgx_get_dose_response_points` to retrieve raw PSet-backed concentration and
  viability points by dataset, drug, and sample.
- `pgx_plot_dose_response` to create a PNG curve plot for a sample-drug pair.
- `pgx_compute_dose_response_metrics` for user-supplied concentration and
  viability vectors.

Mention that fitted metrics depend on dose range, curve quality, and whether
the response crosses the target effect level. Avoid overinterpreting IC50 when
the curve does not cross 50% viability in the tested range.

## Combination And Synergy Analysis

When users ask about drug combinations, Bliss, HSA, ZIP, Loewe, or synergy,
first distinguish:

- raw monotherapy viability
- raw combination viability
- fitted monotherapy dose-response profiles
- fitted combination response surfaces
- summary synergy scores

Use:

- `pgx_list_treatment_response_tables` to find available combo endpoints and
  numeric score columns.
- `pgx_rank_synergy_combinations` to rank high-confidence drug combinations
  from combo treatmentResponse tables.
- `pgx_get_combo_response_records` to inspect selected combination rows.
- `pgx_plot_synergy_heatmap` to create a PNG dose-matrix heatmap for a
  selected or top-ranked combo/sample.
- `pgx_combo_biomarker_association` to associate molecular features with a
  selected combo response or synergy phenotype.
- `pgx_compare_mono_combo_biomarkers` to compare selected feature associations
  for a combo phenotype against each monotherapy arm.
- `pgx_detect_assay_mode` to check whether a dataset appears monotherapy,
  combination, mixed, or perturbation-only.
- `pgx_compute_synergy_reference` only when monotherapy viability vectors are
  supplied or extracted.

Do not imply that real combination data were retrieved unless matched mono and
combination data were actually present. Standard bundled small datasets should
not be treated as combination datasets unless combination assays are explicitly
detected.

For real combination PharmacoSets such as NCI-ALMANAC or GDSC-square, raw
combination dose matrices are required before synergy heatmaps or synergy
scores can be interpreted.

## PSet Curation

When users want help creating or validating a PharmacoSet, use a staged
workflow.

Use:

- `pgx_pset_curation_questions` to ask what tables and columns they have.
- `pgx_validate_pset_inputs` to check whether minimum sample, treatment,
  response, dose, and metadata columns are present.

Do not attempt automatic PharmacoSet construction unless the user has supplied
validated tables and explicitly asks for construction.

## Biomarker And Feature-Response Association

When users ask whether a biomarker, molecular feature, or metadata feature is
associated with drug response, identify:

- dataset or PharmacoSet
- treatment or drug
- sensitivity metric
- feature type: expression, mutation, CNV, pathway score, metadata, etc.
- whether the feature is continuous or categorical

Use:

- `pgx_get_molecular_profile` to inspect molecular values for selected features
  and samples.
- `pgx_suggest_biomarker_candidates` only to seed expected positive-control
  features such as SLC29A1/hENT1 for gemcitabine or PDE3A for anagrelide.
- `pgx_association_test` for exploratory feature-response association.
- `pgx_multi_pset_association_test` to run independent per-PSet associations
  and summarize cross-PSet support without pooling raw data.

Default to Spearman for continuous molecular features versus continuous
sensitivity metrics. Use Wilcoxon/group comparison for binary or categorical
features. Report sample size, effect size, p-value, and adjusted p-value when
available.

Do not claim molecular biomarker associations unless molecular profile data
were actually used. For bundled small datasets, avoid strong biomarker claims;
describe results as exploratory fixtures.

For multi-PSet biomarker questions, do not merge raw response or molecular
matrices across PSets. Run each PSet separately and summarize effect direction,
sample size, p-value/FDR, and consistency.

## Cross-PharmacoSet Analysis

Use:

- `pgx_find_pset_overlap` to find shared samples, treatments, and
  sample-treatment pairs.
- `pgx_compare_pset_response` to compare response values for the same drug and
  metric across datasets.
- `pgx_multi_pset_association_test` for cross-PSet biomarker support using
  per-PSet association followed by consensus summarization.

Check identifier consistency before interpreting cross-PSet results. Matching
names do not always guarantee the same biological sample or compound.

For cross-PSet biomarker support, treat the analysis as reproducibility
evidence, not a single association test. Start with a drug-anchored question and
confirm shared samples, treatments, metrics, and molecular features.

## Limits

- Do not present outputs as clinical treatment recommendations.
- Do not invent data, fitted curves, synergy scores, or biomarkers.
- Do not use arbitrary R execution when a typed MCP tool can answer the
  question.
- State clearly when an analysis is hypothetical, toy-data-backed, or
  PSet-backed.
- If MCP tools do not support the requested workflow, explain what structured
  data or tool support is missing.
