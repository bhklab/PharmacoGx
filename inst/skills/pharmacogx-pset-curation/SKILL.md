---
name: pharmacogx-pset-curation
description: Build, repair, or review PharmacoSet curation pipelines with PharmacoGx and CoreGx, preserving source identifiers, assay semantics, provenance, and scientific QC. Use for source-to-PSet work and PSet structural validation.
---

# PharmacoSet curation

Build a reproducible source-to-PSet workflow, using direct R and the actual
installed PharmacoGx/CoreGx APIs. MCP is optional. Preserve existing inputs and
working curation environments; make migrations explicit and save new versions.

Read [PSet structure](references/pset-structure.md) when designing or repairing
the object, and [scientific QC](references/scientific-qc.md) before deciding
that the curation is complete.

## Source and environment contract

- Record the requested package source/commit, installed versions, R version and
  dependency environment. Check exported constructors and method signatures;
  a newer version number alone does not establish compatibility. Use an exact
  approved revision when development functionality is required.
- Inventory source releases, records, filenames, checksums, formats, licenses
  and measurement definitions. Prefer processed source tables when they retain
  the information needed for the task. Missing processed data is a blocker to
  explain, not permission to silently download bulk raw CEL or sequencing data.
- Preserve source-native keys and original values; normalize via explicit,
  auditable mappings. Display-name equality is not chemical or sample identity.
  Keep unresolved and ambiguous mappings visible instead of guessing.
- Make acquisition, transformation, construction and validation reproducible
  from declared inputs and parameters. Use the project's existing workflow and
  environment tools; Snakemake, Pixi, directory layout, TSV versus CSV, and
  annotation services are conventions rather than PSet validity requirements.

For existing ORCESTRA releases, consult the companion analysis skill's
[dataset handling](../pharmacogx-agentic-analysis/references/dataset-handling.md)
when available. Preserve biological model type and disease context in metadata.
For PRISM/CCLE augmentation, retain response and molecular sources separately,
validate mappings and write a new object with unchanged PRISM response values.

## Essential data invariants

1. Define sample, treatment, feature and assay-level keys and their cardinality.
   Check uniqueness at the intended grain and referential integrity across
   metadata, molecular assays and response tables. Preserve one-to-many assay
   replicates; never equate assay column IDs with biological sample IDs by position.
2. Declare what each assay measures: raw observations, normalized measurements,
   fitted parameters, published endpoints, recomputed endpoints or synergy.
   Preserve concentration units/log base, viability scale, exposure time,
   plate/batch/replicate context, treatment arms and dose pairing.
3. Keep published and recomputed values distinct. Missing, failed fit, censored,
   below detection and true zero are different states. Do not clip, impute or
   drop outliers without a declared rule and an exclusion audit.
4. For molecular profiles, retain source feature identifiers, mapping release,
   assay technology and transformation. Ensembl/Cellosaurus/PubChem mappings are
   useful when appropriate, not universal prerequisites. Do not discard probes,
   transcripts or non-gene features solely to impose gene-level naming.
5. Use public constructors/accessors; inspect legacy versus TRE/MAE storage
   before choosing a path. Validate keys before construction, then run object
   validity checks, accessor-based round trips, source-value comparisons and
   representative analyses. Structural validity alone is insufficient.

## Deliverable and completion

Produce a portable serialized PSet plus the reproducible construction scripts,
input/output manifests, environment provenance and integrated scientific QC.
Keep bulky audit tables and machine-specific paths outside object annotations;
retain compact scientific provenance and mappings in appropriate metadata.
Write new output files, verify reload and content, then promote only the intended
deliverable. Do not mutate published inputs during validation.

QC must explore the actual data, demonstrate representative use, and compare
published versus recomputed results where comparable source data exists.
Report unavailable comparisons explicitly. Use one coherent dataset QC report;
do not create a separate smoke-QC report. Lightweight checks belong in tests.
Human benchmark scripts, expected results and evaluator notes stay outside
agent inputs; generic curation instructions must not encode their answers.
