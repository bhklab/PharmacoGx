---
name: pharmacogx-agentic-analysis
description: Plan and execute reproducible pharmacogenomic analyses with PharmacoGx in direct R, including response curves, biomarkers, cross-PSet comparisons, combinations and study reproduction. Use optional MCP only when that execution route is chosen.
---

# PharmacoGx agentic analysis

Use this Skill as a portable scientific workflow layer for an external host
such as Biomni or Claude. It is not a standalone agent.

## Choose the execution route

Prefer direct R; select one route before analyzing data and report it:

- **Direct R:** read [references/direct-r.md](references/direct-r.md) for
  version verification, safe execution and package-specific decisions. In
  Biomni web also read [references/biomni-web.md](references/biomni-web.md).
- **MCP:** when explicitly chosen and the host exposes typed `pgx_*` tools,
  read [references/mcp-tools.md](references/mcp-tools.md).
- **Planning only:** if neither route is available, return a reproducible plan
  and the missing capabilities. Do not imply that an analysis was executed.

Do not silently switch routes. A Skill-only direct-R run is not a Skill+MCP
plugin run.

For response summaries/curves, molecular associations, cross-PSet comparisons,
combination analysis or paper reproduction, read
[analysis workflows](references/analysis-workflows.md). A compact
[development example](scripts/direct-r-examples.R) uses only bundled/synthetic
data; it is an API check, not scientific validation or a held-out benchmark.

For source-to-PSet construction, use the companion `pharmacogx-pset-curation`
Skill when available. For remote measurement/catalog discovery, use
`orcestra-discovery`. Neither is required to analyze an already verified local
PSet. Do not expose evaluator scripts or expected results to benchmark agents.

## Dataset-aware handling

Before selecting cohorts, molecular profiles or combination summaries, read
[dataset handling](references/dataset-handling.md). It covers biological model
types, tissue-specific PSets, TRE score selection and known dataset differences.
For PRISM with missing molecular measurements, read
[PRISM and CCLE integration](references/prism-ccle.md) before proposing a
validated pairing. External molecular augmentation must fit the user's scope
and any frozen benchmark protocol; it is a distinct data condition.

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

Identify whether samples are established cell lines, primary patient specimens,
patient-derived cultures or xenograft-derived cultures. Preserve that distinction
in summaries and use consistent model-type colors and a legend in cross-PSet
plots. Consider relevant disease-specific datasets even when tissue labels are
missing or use equivalent terminology; verify source provenance and retain
contradictory sample annotations for review. A frozen cohort definition takes
precedence over a generic dataset recommendation.

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
