# PharmacoGx Agentic Integration Demo

This directory contains an optional Model Context Protocol (MCP) demo for
calling a small, curated PharmacoGx analysis surface from agentic systems.

The demo is intentionally narrow. It exposes JSON-serializable wrappers around
bundled toy datasets and dose-response utilities rather than raw S4 objects or
arbitrary R evaluation.

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
- "Compute AUC, AAC, IC50, and AC50 for this dose-response vector."
- "Compute Bliss and HSA references for these two monotherapy viabilities."
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

- Prefer bundled demo datasets for reproducible agent demos.
- Call `pgx_list_available_covariates()` before metadata lookup or filtering,
  because annotation fields differ between PharmacoSets.
- Use `pgx_list_available_psets()` to discover downloadable external datasets.
- Do not call `pgx_download_pset()` unless the user explicitly asks for a large
  external dataset workflow and the tool call sets `confirm_download = TRUE`.
- Treat outputs as preclinical pharmacogenomic analysis, not clinical treatment
  advice.
- Keep returned tables small enough for agent context windows.
