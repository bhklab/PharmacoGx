# PharmacoGx Agentic Integration Demo

This directory contains an optional Model Context Protocol (MCP) demo for
calling a small, curated PharmacoGx analysis surface from agentic systems.

The demo exposes JSON-serializable wrappers around bundled toy datasets, local
`.qs`/`.rds` PharmacoSet files, and dose-response utilities rather than raw S4
objects or arbitrary R evaluation.

## Local PSet Manifest

The MCP tools also look for `inst/mcp/local-pset-manifest.csv`, or the path in
`PHARMACOGX_LOCAL_PSET_MANIFEST`, to resolve short dataset aliases to local
`.qs`/`.rds` PharmacoSet files. In the PGX3 development workspace, the manifest
includes:

- `nci_almanac`: canonical NCI-ALMANAC combo PSet at
  `/Users/michael/Projects/BHKLab/pgx3/nci_almanac_pset_5.rds`.
- `gdsc2_matrix`: GDSC2 Matrix combo PSet.
- `gdsc2_anchor`: GDSC2 Anchor combo PSet.
- `prism`, `hmcl`, and `nci60_legacy`: additional local demo/reference PSets.

Call `pgx_list_local_psets()` before real-data demos to confirm which aliases
exist on the current machine. Explicit local file paths still work when a PSet
is not listed in the manifest.

The local aliases are read-only analysis inputs. `nci_almanac` is large and can
take about a minute to load from cold storage, so use one long-lived R/MCP
session for multi-step combo demos. Some legacy aliases, such as `prism` and
`nci60_legacy`, use older serialized PharmacoSet slots; the MCP includes
read-only fallbacks for metadata, molecular profile, and sensitivity-profile
access where those legacy tables are present.

## Requirements

Install the optional agent integration packages in the R environment used by
your MCP client:

```r
if (!requireNamespace("pak", quietly = TRUE)) {
  install.packages("pak")
}
pak::pkg_install(c("mcptools", "ellmer", "btw"))
```

`btw` is optional but useful when an MCP client also needs general R session
inspection tools.

## Start The MCP Server

From a source checkout:

```sh
Rscript inst/mcp/run-pharmacogx-mcp.R
```

When the runner is started from a source checkout, it uses the sibling
`PharmacoGx-tools.R` file before any installed-package copy. Restart the MCP
server after editing the tool file; an already-running stdio server keeps the
version it sourced at startup. The runner also loads the PharmacoGx namespace
from the checkout with `pkgload` by default, avoiding mismatches between the
branch and an older installed package. Set
`PHARMACOGX_MCP_USE_SOURCE_PACKAGE=false` to use the installed package instead.

Set `PHARMACOGX_MCP_LOG_FILE` to capture tool start, success, and R error
details without writing diagnostics to the stdio protocol stream. The VS Code
example in `.vscode/mcp.json` writes to `.vscode/pharmacogx-mcp.log`.
`pgx_session_info()` reports the active tool file, its modification time, the
server start time, process ID, and log path for troubleshooting stale servers.

From an installed package:

```sh
Rscript -e "mcptools::mcp_server(tools = system.file('mcp', 'PharmacoGx-tools.R', package = 'PharmacoGx'))"
```

Example client configuration:

```json
{
  "mcpServers": {
    "pharmacogx": {
      "command": "Rscript",
      "args": ["/absolute/path/to/PharmacoGx/inst/mcp/run-pharmacogx-mcp.R"]
    }
  }
}
```

To let the MCP server route tool calls into a live R session, run this inside
that R session:

```r
library(PharmacoGx)
mcptools::mcp_session()
```

If no live session is registered, the tools run in the MCP server process.

## Demo Prompts

- "List the PharmacoGx example datasets and what each supports."
- "Using GDSCsmall, list available sensitivity measures."
- "Rank the top Doxorubicin responders in GDSCsmall."
- "Why might the top Doxorubicin responders be similar? Check sample metadata
  before answering."
- "How robust is the Doxorubicin responder ranking across available sensitivity
  measures?"
- "Find GDSCsmall samples from a urogenital or breast lineage and summarize
  their metadata."
- "Retrieve and plot the GDSCsmall Doxorubicin dose-response curve for 22RV1."
- "Compute AUC, AAC, IC50, and AC50 for this dose-response vector."
- "Compute Bliss and HSA references for these two monotherapy viabilities."
- "Ask me the information needed to curate a sensitivity PharmacoSet from my
  tables, then validate the columns I provide."
- "Test whether ALK mutation or TSPAN6 expression is associated with
  Doxorubicin response in the toy data, and explain the limitations."
- "Which samples, treatments, and sample-treatment pairs overlap between
  GDSCsmall and CCLEsmall?"
- "List local PSet aliases, then use `nci_almanac` to list treatmentResponse
  tables and score candidates."
- "Using `nci_almanac`, rank high-confidence synergistic combinations by
  `ZIP_delta`, then plot a heatmap for the top pair."
- "Using `nci_almanac`, compare selected RNA biomarker associations for a drug
  pair against each monotherapy arm."
- "Using this local GDSC-square Matrix `.qs` file path, list treatmentResponse
  tables and rank high-confidence synergistic drug combinations."
- "Find candidate biomarkers for Gemcitabine, then test SLC29A1/hENT1 support
  across available PSets without pooling raw data."
- "For this GDSC-square drug combo, test whether mutation or CNV features are
  associated with the combo synergy score."
- "List remote PharmacoSets available for download, then ask me before
  downloading one."

## Biomni

Biomni can import external MCP servers with `agent.add_mcp(config_path = ...)`.
Use `biomni-mcp-config.yaml` as a template and replace the script path with the
absolute path to `run-pharmacogx-mcp.R`.

## ToolUniverse

ToolUniverse's MCP auto-loader expects an HTTP MCP server. Start this demo in
HTTP mode first:

```sh
PHARMACOGX_MCP_TRANSPORT=http MCPTOOLS_PORT=8080 \
  Rscript inst/mcp/run-pharmacogx-mcp.R
```

Then load `tooluniverse-mcp-tools-config.json` from ToolUniverse.

## Agent Skill

The companion skill is at:

```text
inst/skills/pharmacogx-agentic-analysis/SKILL.md
```

Install or copy that skill into the skill directory supported by your agent
harness. The skill is an instruction/routing layer only; MCP remains the
execution layer.

## Guardrails

- Prefer bundled demo datasets for reproducible agent demos; use explicit local
  `.qs` or `.rds` file paths for real PSet-backed workflows.
- Clearly distinguish toy-data demonstrations from real downloaded PSet-backed
  analyses.
- Call `pgx_list_available_covariates()` before metadata lookup or filtering,
  because annotation fields differ between PharmacoSets.
- Use `pgx_detect_assay_mode()` before any drug-combination or synergy
  workflow.
- Use `pgx_list_treatment_response_tables()` before combo ranking or combo
  biomarker workflows, because combo endpoints may live in `profiles`,
  `synergy`, or `raw` treatmentResponse tables.
- For multi-PSet biomarker workflows, run per-PSet associations and summarize
  consistency; do not pool raw response or molecular matrices across PSets.
- Use `pgx_list_available_psets()` to discover downloadable external datasets.
- Do not call `pgx_download_pset()` unless the user explicitly asks for a large
  external dataset workflow and the tool call sets `confirm_download = TRUE`.
- Treat outputs as preclinical pharmacogenomic analysis, not clinical treatment
  advice.
- Keep returned tables small enough for agent context windows.
