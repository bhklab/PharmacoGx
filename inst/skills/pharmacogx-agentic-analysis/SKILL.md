---
name: pharmacogx-agentic-analysis
description: Plan and execute reproducible preclinical pharmacogenomic analyses with PharmacoGx from a host agent, using typed MCP tools when available or a constrained direct-R workflow in environments such as Biomni web.
---

# PharmacoGx agentic analysis

Use this Skill as a portable scientific workflow layer for an external host
such as Biomni or Claude. It is not a standalone agent.

## Choose the execution route

Select one route before analyzing data and report it in the result:

- **MCP:** if the host exposes `pgx_session_info` and the other typed
  `pgx_*` tools, read [references/mcp-tools.md](references/mcp-tools.md).
- **Direct R:** if the host can execute R with PharmacoGx but has no
  PharmacoGx MCP tools, read
  [references/biomni-web.md](references/biomni-web.md). This is the immediate
  Biomni web route. For the current pre-merge demo, that reference pins
  PharmacoGx to the GitHub `feat/mcp` branch and requires source verification
  before analysis.
- **Planning only:** if neither route is available, return a reproducible plan
  and the missing capabilities. Do not imply that an analysis was executed.

Do not silently switch routes. A Skill-only direct-R run is not a Skill+MCP
plugin run.

## Workflow contract

1. Identify the requested treatment, sample or disease context, response
   endpoint, molecular profile, feature type, and candidate features. Resolve
   synonyms only as needed and preserve stored identifiers.
2. Check whether candidate PharmacoSets jointly contain every required element
   before analysis. Distinguish compatible, incompatible, unchecked, and
   load-failed datasets. A catalog entry or named molecular profile is not
   proof that usable measurements exist.
3. Use bundled fixtures only for tests and lightweight examples. Label toy and
   real PSet-backed results explicitly.
4. Before a remote download, show a filled manifest with the exact file,
   source record, version, size, expected checksum, persistent storage
   location, and compatibility rationale. Ask for approval after showing that
   manifest. A dataset choice made during planning is not download approval.
5. If compatibility cannot be checked without downloading, say so and propose
   the smallest staged download. Do not download a batch of unchecked PSets
   merely because they appear in a distribution catalog.
6. Inspect actual sample, treatment, sensitivity, and every relevant molecular
   profile before choosing fields. Compare usable cohort and complete-case
   coverage using a pre-specified rule; do not select the first matching
   profile.
7. Freeze the cohort, endpoint, molecular profile, identifier mappings,
   analysis method, inference settings, and multiplicity family before
   inspecting the association result. Report initial matches, exclusions,
   missingness, and final complete-case counts.
8. Keep outputs bounded and preserve a complete analysis script plus compact
   machine-readable results. Follow
   [references/reporting-contract.md](references/reporting-contract.md).

## Scientific rules

### Drug response

State the metric definition and direction before ranking or associating it:

- Lower AUC usually indicates greater sensitivity.
- Lower IC50 usually indicates greater sensitivity when it is identified
  within the tested concentration range.
- Higher AAC-like or activity-area values usually indicate greater
  sensitivity.

Compare response metrics before claiming that responder rankings are robust.
Treat small cohorts and singleton metadata groups as exploratory.

### Sample context

Use the exact available annotation field and stored value. Tissue, primary
site, lineage, and disease or histology are not interchangeable. For example,
do not label all pancreas-lineage models as pancreatic ductal adenocarcinoma
without supporting disease or histology metadata.

Make the primary cohort match the requested disease context when authoritative
disease or histology annotations permit it. If only a broader tissue or lineage
cohort can be constructed, rename the analysis and every output to that broader
context. Analyze a broader cohort only as a separately labeled sensitivity
analysis.

Metadata-only associations are not molecular biomarkers.

### Biomarkers

For a continuous molecular feature and continuous response, default to a
rank-based association unless the user or a frozen protocol specifies another
model. For categorical features, use a suitable group comparison or regression
model. Report sample size, effect size, direction, uncertainty, raw p-value,
and adjusted p-value when applicable.

Define one primary endpoint and method before execution. Define the family of
confirmatory tests across datasets, endpoints, and features, and either adjust
within that family or give a specific pre-specified reason for not doing so.
Label secondary endpoints and alternative methods as sensitivity or
exploratory analyses. Do not infer that adjustment is unnecessary merely
because only one molecular feature was requested.

Suggestions from prior knowledge are hypotheses, not evidence that a feature
exists in the selected PSet or is associated with the treatment. Do not claim
validation from weak, directionally conflicting, or single-PSet evidence.

For multiple PSets, analyze each independently and summarize signed effects,
uncertainty, sample counts, endpoint definitions, and feature mappings. Do not
pool raw response or molecular values across PSets without a pre-specified
integration model and justification. Do not assign discovery and validation
roles after seeing the results.

### Combination response

First determine whether the data contain monotherapy responses, combination
responses, fitted surfaces, or summary synergy scores. Do not calculate or
interpret Bliss, HSA, ZIP, or Loewe results unless the required inputs are
present and their meaning is documented. Distinguish a combination-specific
association from an association driven by either monotherapy arm.

### PSet curation

Before constructing a PharmacoSet, identify and validate the sample,
treatment, response, dose, molecular, and annotation tables. Do not construct
an object until the user has supplied the required tables and requested that
operation.

## Evidence boundary

- Separate executed PSet evidence, attached capability metadata, and external
  biological literature.
- Treat a failed feasibility check or a well-supported negative result as a
  valid outcome.
- Do not present preclinical findings as clinical treatment recommendations.
- Do not invent measurements, identifiers, fitted curves, synergy scores, or
  biomarkers.
- Explain unavailable operations and required inputs instead of hiding them
  behind an unrecorded fallback.
