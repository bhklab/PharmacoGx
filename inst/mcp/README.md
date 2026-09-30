# PharmacoGx Agent Plugin

This directory contains the typed execution layer for a portable PharmacoGx
plugin. External hosts load the shared PharmacoGx Skill for workflow policy.
Hosts with MCP support can call this server for bounded execution; hosts such
as Biomni web can initially use the same Skill with a declared direct-R route.
PharmacoGx does not ship a standalone agent or its own model-level planner.

The demo exposes JSON-serializable wrappers around bundled toy datasets, local
`.qs`/`.rds` PharmacoSet files, and dose-response utilities rather than raw S4
objects or arbitrary R evaluation.

## Local PSet Discovery

No developer machine, shared lab mount, or pre-existing data collection is
required. Each user runs R and this MCP server on their own computer or an
approved compute host. The model host calls the server; it does not execute R
merely by loading the Skill. A hosted, authenticated remote MCP service is a
separate deployment, not a prerequisite for this local stdio MVP.

### Downloading from the current ORCESTRA service

1. Call `pgx_list_available_psets(canonical = FALSE)` to obtain the current
   ORCESTRA catalogue and exact PSet names.
2. Ask the user to approve the dataset, download location and resource use.
3. Call `pgx_download_pset(name = "<exact PSet name>", confirm_download = TRUE)`.
   This delegates to `PharmacoGx::downloadPSet()`. By default it saves under
   `file.path(tools::R_user_dir("PharmacoGx", "data"), "psets")`, not a temporary
   directory. This directory is automatically scanned on subsequent sessions.
4. Pass the returned `file_path` to `pgx_check_analysis_requirements()` and
   analysis tools. A custom download directory need not be in the discovery
   roots; use the file path directly or configure `PHARMACOGX_LOCAL_PSET_DIRS`.

The existing ORCESTRA catalogue supplies the download URL. Downloads may still
use Zenodo and be slow; increase the timeout when appropriate. Neither entity
coverage discovery (genes, cell lines, drugs) nor GCP-backed downloads is assumed
available. Both remain future ORCESTRA collaboration work. Catalogue membership
alone does not establish that a PSet supports a particular analysis. No large
download is made automatically to answer a feasibility question.

`pgx_list_local_psets()` combines two local discovery mechanisms:

1. An optional manifest from `PHARMACOGX_LOCAL_PSET_MANIFEST` or
   `inst/mcp/local-pset-manifest.csv`.
2. Recursive discovery of likely PSet `.rds` and `.qs` files under the
   directories in `PHARMACOGX_LOCAL_PSET_DIRS`, plus the source checkout's
   `local_data` directory and the PharmacoGx user data directory.

Multiple discovery directories use the platform path separator. Manifest paths
may be relative to the manifest file, so shared configuration does not require
machine-specific absolute paths. Explicit local file paths also remain valid.

Auto-discovered names remove a leading `PSet_` and use lowercase snake case.
Known local names such as `NCI_ALMANAC_2017.rds` and `PSet_PRISM.rds` resolve to
`nci_almanac` and `prism`. Use `pgx_check_analysis_requirements()` before
loading several PSets or attempting an association.

Local requirements checks run in isolated R workers and retain only compact
capability snapshots in the long-lived MCP process. A large or malformed PSet
therefore cannot consume the server's memory and cause later datasets to be
reported as unavailable. Snapshots are invalidated when the PSet or MCP tool
file changes and are cached under `PHARMACOGX_MCP_CAPABILITY_CACHE_DIR`, or the
platform PharmacoGx cache directory by default. Set
`PHARMACOGX_MCP_PROBE_TIMEOUT` to change the default 300-second per-PSet probe
limit.

The requirements check resolves exact identifiers across PSet annotations,
including treatment names, PubChem CIDs, gene symbols, versioned Ensembl IDs,
and dataset-specific row identifiers. When no PSet meets every requirement, it
returns the closest successfully checked candidates and separates them from
files that could not be loaded. The remote `availablePSets()` table is not a
coverage index and is not used as an automatic fallback.

Directories named `PSets` are treated as curated PSet collections. In broader
project trees, discovery requires `pset` in the filename or a recognized PSet
filename, which prevents unrelated serialized analysis tables from entering
the index.

## Requirements

Install the optional agent integration packages in the R environment used by
your MCP client:

```r
if (!requireNamespace("pak", quietly = TRUE)) {
  install.packages("pak")
}
pak::pkg_install(c("bhklab/PharmacoGx@feat/mcp", "mcptools", "ellmer"))
```

Until this branch is merged, install `feat/mcp`, not the release package. Record
the installed `RemoteSha` for reproducibility. For checkout development also
install `pkgload`; transport smoke testing additionally requires `processx`.
Local unpushed changes are not included in a GitHub installation.

General R session tools are disabled by default to keep VS Code focused on the
curated PharmacoGx surface. Set `PHARMACOGX_MCP_SESSION_TOOLS=true` only when a
live general-purpose R session is intentionally needed.

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
Rscript -e "source(system.file('mcp', 'run-pharmacogx-mcp.R', package = 'PharmacoGx'), chdir = TRUE)"
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

## Development Host Setup

The workspace includes `.vscode/mcp.json` for VS Code. Use the client example
above for other MCP-compatible development hosts. These are adapters for exercising
the same plugin surface that Claude, Biomni, and other compatible hosts load.
The VS Code configuration:

- starts the source runner with `${workspaceFolder}` paths
- scans `${workspaceFolder}/local_data` and the current user's PSet data directory
- returns one JSON content block per tool call
- disables unrelated R-session tools
- restarts the MCP server when an `inst/mcp/*.R` file changes

After opening the checkout in VS Code:

1. Run **MCP: List Servers**, select `pgx`, and choose **Restart**.
2. If tools were previously cached, run **MCP: Reset Cached Tools**.
3. Ask a pharmacogenomic question through the host's ordinary agent. For
   biomarker questions, the
   generic skill can load automatically, or invoke it explicitly with
   `/biomarker-validation` followed by the question.
4. Confirm that the first call is `pgx_session_info` and that a requirements
   preflight occurs before analysis. Treat terminal or generic-R fallback as an
   evaluation failure when the question requires plugin-backed evidence.

The reusable instruction layers are:

```text
.github/skills/biomarker-validation/SKILL.md
inst/skills/pharmacogx-agentic-analysis/SKILL.md
```

Run the transport-level smoke test after changing the registry or runner:

```sh
Rscript inst/mcp/smoke-test-pharmacogx-mcp.R
```

To intentionally expose general R-session routing, set
`PHARMACOGX_MCP_SESSION_TOOLS=true`, restart the server, and run this inside the
target R session:

```r
library(PharmacoGx)
mcptools::mcp_session()
```

The normal VS Code configuration keeps this disabled. Analysis tools share one
persistent MCP process, while local capability discovery uses isolated workers
and stores only lightweight snapshots in the server.

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
- "For a biomarker, treatment, and disease context I provide, identify which
  PSets contain all required response, sample, and molecular data before
  attempting an association."
- "Validate a user-provided molecular feature against treatment response in
  each compatible PSet independently, then compare effect direction and
  uncertainty without pooling raw data."
- "For this GDSC-square drug combo, test whether mutation or CNV features are
  associated with the combo synergy score."
- "Use `availablePSets()` through the remote catalog tool to list downloadable
  PharmacoSets, then ask me before calling the guarded `downloadPSet()` tool."

## Biomni

Biomni has two distinct integration paths:

1. The open-source Python client can start a local stdio MCP server with
   `agent.add_mcp(config_path = ...)`. Use `biomni-mcp-config.yaml` as a
   template and replace the script path with the absolute path to
   `run-pharmacogx-mcp.R`.
2. Biomni web can attach the portable Skill and execute PharmacoGx through R.
   Its custom MCP connector requires an organization administrator to register
   a remotely reachable HTTPS MCP URL and configure OAuth or per-user token
   authentication. The local stdio template is not a web-portal connector.

For the immediate Biomni web route, attach the complete
`inst/skills/pharmacogx-agentic-analysis` folder and follow
`references/biomni-web.md`. Record the evaluation condition as `skill`, not
`plugin`. During pre-merge testing, that route installs exactly
`bhklab/PharmacoGx@feat/mcp` with `pak` and verifies the installed Git reference
and commit before analysis. It must not fall back to Bioconductor or another
GitHub branch. A future hosted, authenticated server can be evaluated as
`plugin`.

## ToolUniverse

ToolUniverse's MCP auto-loader expects an HTTP MCP server. Start this demo in
HTTP mode first:

```sh
PHARMACOGX_MCP_TRANSPORT=http MCPTOOLS_PORT=8080 \
  Rscript inst/mcp/run-pharmacogx-mcp.R
```

Then load `tooluniverse-mcp-tools-config.json` from ToolUniverse.

## Plugin Skill

The companion skill is at:

```text
inst/skills/pharmacogx-agentic-analysis/SKILL.md
```

Install or copy the complete skill folder into the skill area supported by the
host platform so its conditional references remain available. The Skill is the
shared workflow layer. MCP is the preferred typed execution route when it is
available, while direct R is an explicit constrained route for hosts without
the connector. Host-specific configuration should remain a thin adapter.

## Evaluation

The host-neutral development suite lives under `inst/agentic-eval/`. Validate
the case manifest and initialize a run record with:

```sh
Rscript inst/agentic-eval/scripts/validate-cases.R
Rscript inst/agentic-eval/scripts/init-run.R \
  --host biomni-web \
  --condition skill \
  --model MODEL_ID \
  --plugin-commit GIT_COMMIT \
  --output /path/to/run.json
```

Claude, Biomni, and future hosts should receive identical prompts for matched
comparisons. Host adapters write the raw response, execution trace, artifacts,
timing, and adjudication into the initialized record. Development cases may be
used to improve the plugin; locked benchmark cases must not be used for tuning.

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
