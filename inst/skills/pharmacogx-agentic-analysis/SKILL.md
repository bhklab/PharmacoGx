---
name: pharmacogx-agentic-analysis
description: Use PharmacoGx for preclinical pharmacogenomic analysis with MCP-backed R workflows.
---

## Purpose

Use this skill when a user asks an agent to perform or explain preclinical
pharmacogenomic analysis with PharmacoGx.

Prefer the PharmacoGx MCP tools when they are available. They provide a small,
typed analysis surface for agent workflows:

- `pgx_list_example_datasets`
- `pgx_list_entities`
- `pgx_list_available_covariates`
- `pgx_get_sample_metadata`
- `pgx_get_treatment_metadata`
- `pgx_filter_samples`
- `pgx_compare_metrics`
- `pgx_summarize_sensitivity`
- `pgx_top_responders`
- `pgx_compute_dose_response_metrics`
- `pgx_compute_synergy_reference`
- `pgx_list_available_psets`
- `pgx_download_pset`
- `pgx_session_info`

## Workflow

1. Start by listing datasets, entities, or covariates if the user's requested
   drug, cell line, metadata field, or sensitivity measure is ambiguous.
2. Use bundled toy datasets (`GDSCsmall`, `CCLEsmall`, `CMAPsmall`) by default.
   Use `pgx_list_available_psets` for external dataset discovery only when the
   user asks for larger downloadable PharmacoSets.
3. Call `pgx_list_available_covariates` before metadata lookup or filtering,
   because annotation field names differ between PharmacoSets.
4. Keep tool outputs small. Request ranked or filtered results instead of full
   sensitivity matrices.
5. Explain the metric direction when ranking response:
   - Lower AUC or IC50 usually indicates stronger sensitivity.
   - Higher AAC-like response area usually indicates stronger sensitivity.
6. Use `pgx_compare_metrics` before claiming a responder pattern is robust.
7. Include package/session metadata for reproducibility when summarizing a run.

## Limits

- Do not present PharmacoGx outputs as clinical treatment recommendations.
- Do not call `pgx_download_pset` unless the user explicitly asks for that
  workflow and accepts the runtime, storage, and network cost. The tool must be
  called with `confirm_download = TRUE`.
- Do not use arbitrary R execution when a typed PharmacoGx MCP tool can answer
  the question.
- If a requested analysis requires unsupported inputs, state what structured
  data is missing rather than inventing results.

## Setup Reference

The MCP demo assets are in `inst/mcp/`:

- `PharmacoGx-tools.R` defines the MCP tools.
- `run-pharmacogx-mcp.R` starts the server.
- `biomni-mcp-config.yaml` is a Biomni template.
- `tooluniverse-mcp-tools-config.json` is a ToolUniverse HTTP auto-loader
  template.
