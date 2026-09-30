---
name: biomarker-validation
description: Validate whether a molecular biomarker is associated with drug response in preclinical pharmacogenomic datasets using PharmacoGx. Use for biomarker-treatment questions involving expression, mutation, copy number, or other molecular features across one or more PharmacoSets, including feasibility checks, cohort definition, association testing, cross-dataset comparison, and evidence limitations.
---

# Biomarker validation

Use this workflow for the biomarker, treatment, disease context, and response
endpoint supplied by the user. Do not substitute entities from examples or
assume a dataset contains the requested data.

1. Choose and disclose the execution route. If PharmacoGx MCP tools are
   available, call `pgx_session_info` once and confirm the expected source
   checkout and tool registry. If only direct R is available, record the R and
   PharmacoGx versions and follow the direct-R route in
   `inst/skills/pharmacogx-agentic-analysis/references/biomni-web.md`.
2. Resolve only necessary synonyms, then check all requested requirements
   together before analysis. With MCP, call `pgx_check_analysis_requirements`
   once. With direct R, use a versioned capability snapshot when available or
   clearly describe the smallest staged inspection required. Do not invent an
   AUC or AAC endpoint when the user did not specify one. Prefer real local or
   verified downloaded PSets over bundled fixtures.
3. Inspect identifier resolution and closest checked candidates before
   interpreting a failed all-requirements check. Treatment names, PubChem CIDs,
   gene symbols, Ensembl IDs, and dataset-specific identifiers may refer to the
   same entity. Preserve the stored identifier and disclose any partial or
   ambiguous mapping.
4. Continue only with datasets that jointly satisfy the requirements. If none
   do, return a capability table identifying each missing element and separate
   datasets that were checked from datasets that could not be loaded. Do not
   jump directly to `availablePSets()`: its catalog does not establish drug,
   disease, endpoint, or molecular coverage. Never download without explicit
   approval.
5. Inspect available covariates before choosing fields. Define the final cohort
   using an explicit tissue or disease field and exact matching. Report the
   field, stored value, exclusions, and resulting sample count.
6. Verify the stored treatment and feature identifiers. Discover the response
   endpoint rather than assuming an AUC-like field exists, and state whether
   larger or smaller values represent greater sensitivity.
7. Run the association only after the cohort, endpoint, molecular profile, and
   identifiers are fixed. Report sample size, effect size, direction, raw and
   adjusted p-values, missingness, and any transformations or summaries.
8. For multiple PSets, analyze each independently and compare direction and
   magnitude. Do not pool values across studies without an explicit
   harmonization method and justification.
9. Separate evidence from interpretation. A failed feasibility check is a
   valid result. Preclinical cell-line evidence does not establish clinical
   validity.

Keep outputs bounded. Do not replace an unavailable typed tool with hidden
terminal or generic-R execution. Direct R is acceptable only when it was the
declared route and the complete script is retained.
