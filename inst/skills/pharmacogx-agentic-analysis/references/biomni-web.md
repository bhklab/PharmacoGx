# Biomni web direct-R execution

Use this route when a Biomni web task can execute R but does not have the
PharmacoGx MCP connector. The Skill supplies workflow policy; PharmacoGx and R
supply execution.

## Temporary package source

Until the MCP work is merged, install PharmacoGx directly from the BHKLab
GitHub `feat/mcp` branch in a fresh R session:

```r
if (!requireNamespace("pak", quietly = TRUE)) {
  install.packages("pak", repos = "https://cloud.r-project.org")
}

pak::pkg_install("bhklab/PharmacoGx@feat/mcp")
```

Do not substitute the Bioconductor release, GitHub `main` or `devel`, or an
unverified preinstalled copy. If this exact installation fails, stop and
report the error instead of silently falling back. If PharmacoGx was already
loaded before installation, restart R before continuing.

Before analysis, verify and record the installed source:

```r
pharmacogx_description <- utils::packageDescription("PharmacoGx")
pharmacogx_description[c(
  "Package",
  "Version",
  "RemoteType",
  "RemoteHost",
  "RemoteUsername",
  "RemoteRepo",
  "RemoteRef",
  "RemoteSha"
)]
```

Require `RemoteUsername` to be `bhklab`, `RemoteRepo` to be `PharmacoGx`, and
`RemoteRef` to be `feat/mcp`. Record `RemoteSha` as the package commit. Treat a
missing or different source field as a failed provenance check.

## Boundaries

- Do not call or describe `pgx_*` tools as available.
- Do not claim that the MCP or ORCESTRA API was used.
- Use only the source-verified PharmacoGx package, attached Skill resources,
  and user-approved public PSet downloads.
- Do not source executable code from an unreviewed remote URL.
- Record every R command or save a complete script. A hidden or irreproducible
  fallback is a failed run.

## Data discovery and downloads

1. Record `R.version.string`, `packageVersion("PharmacoGx")`, the verified
   GitHub source fields and commit, and `sessionInfo()`.
2. Use `availablePSets()` only as a distribution catalog. It does not verify
   whether a PSet contains the requested drug, sample context, response metric,
   or molecular feature.
3. Prefer a versioned capability snapshot attached to the task when it matches
   the exact PSet file version and checksum. Treat entries absent from a partial
   snapshot as unchecked, not unavailable.
4. If no trustworthy capability snapshot exists, state that compatibility is
   unknown. Rank candidates using documented source metadata, propose one
   smallest informative download, audit it, and only then propose another.
   Do not describe catalog entries as compatible or download all unchecked
   candidates together.
5. Before any download, present a completed approval table containing the exact
   filename, PSet version, source record URL, expected bytes, expected checksum,
   persistent location, requirements supported by existing evidence, and
   remaining unchecked requirements. Ask explicitly whether to download those
   listed files. Earlier choices or approval of an analysis plan do not satisfy
   this download gate.
6. Save approved PSets in persistent project or shared storage. Verify actual
   bytes and checksum before loading, reuse the verified file in later steps,
   and retain both expected and observed values in `download_manifest.csv`.

## Freeze the cohort and molecular profile

Before viewing any association result:

1. Resolve the exact requested disease or histology from authoritative sample
   annotations. Use that as the primary cohort when possible. If only tissue or
   lineage is available, use the broader label in the title, manifest, tables,
   figures, and interpretation. Keep any broader cohort as a separately labeled
   sensitivity analysis.
2. Enumerate every molecular profile containing the requested feature. For each
   candidate, record assay name and version, stored feature identifier, matched
   cohort samples, non-missing feature values, matched endpoint values, and
   final complete-case count.
3. Choose the profile using a pre-specified rule, normally the largest valid
   complete-case cohort after assay and identifier checks. Document scientific
   or quality reasons for overriding that rule. Never select a profile merely
   because it is the first match.
4. Save this comparison as `profile_selection.csv` and freeze the chosen cohort,
   profile, endpoint, feature, and identifiers in the analysis manifest.

## Direct PharmacoGx workflow

Load and inspect one large PSet at a time unless the environment has enough
memory for more. Use documented PharmacoGx accessors such as `sampleInfo()`,
`treatmentInfo()`, `sensitivityInfo()`, and `molecularProfiles()` to verify
coverage before analysis. Use `sensitivityMeasures()` or the applicable
sensitivity annotations to confirm the response field.

For a biomarker analysis:

1. Resolve the stored treatment, sample, molecular-profile, and feature
   identifiers.
2. Define the cohort from an observed annotation field and report the selection
   logic.
3. Report counts at capability match, cohort filtering, molecular/response
   intersection, and final complete-case stages.
4. Run `drugSensitivitySig()` only after explicitly fixing `mDataType`,
   `drugs`, `features`, `cells` or `tissues`, `sensitivity.measure`, molecular
   and sensitivity summary statistics, standardization, `modeling.method`,
   `inference.method`, and `returnValues`. Do not rely on positional defaults.
   Extract results by named array dimensions, such as `"estimate"`, `"pvalue"`,
   and `"fdr"`, rather than indexing the internal `.Data` slot.
5. Analyze each PSet separately. Save compact results immediately and release
   the large object before loading the next PSet if memory is constrained.

Define the primary endpoint and method before execution. When more than one
confirmatory dataset, endpoint, or feature is tested, define the multiplicity
family and report adjusted values. Keep secondary endpoints and alternative
methods explicitly exploratory or sensitivity analyses. Report uncertainty for
each effect in a machine-readable table, not only in a plot.

If `drugSensitivitySig()` is intended to reproduce a prior analysis, pass the
prior modeling and inference settings explicitly. If it is intended as a
rank-based check, pass `modeling.method = "spearman"` and a declared
`inference.method`. Inspect the recorded call and arguments before explaining
why it differs from a separate correlation.

Use this pattern after freezing the values on the right-hand side:

```r
dsig_modeling_method <- "spearman"
dsig_inference_method <- "analytic"

dsig <- drugSensitivitySig(
  object = pset,
  mDataType = profile_name,
  drugs = treatment_id,
  features = feature_id,
  cells = primary_cohort_ids,
  sensitivity.measure = endpoint_name,
  molecular.summary.stat = "mean",
  sensitivity.summary.stat = "mean",
  returnValues = c("estimate", "pvalue", "fdr"),
  standardize = "SD",
  modeling.method = dsig_modeling_method,
  inference.method = dsig_inference_method,
  verbose = FALSE
)

dsig_fields <- dimnames(dsig)[[3]]
stopifnot(all(c("estimate", "pvalue") %in% dsig_fields))
dsig_values <- dsig[
  feature_id,
  treatment_id,
  intersect(c("estimate", "se", "pvalue", "fdr"), dsig_fields),
  drop = TRUE
]
```

Replace the two method variables only with values frozen in the analysis plan.
Record them and the complete call in the manifest.

When comparing with an earlier run, preserve its primary endpoint and method
as a matched analysis. Additional methods should be labeled pre-specified
robustness analyses rather than replacements chosen after seeing results. Keep
the prompt, host model, user clarifications, project resources, and approval
decisions identical between evaluation conditions; otherwise label the runs as
illustrative rather than matched.

## Required artifacts

Return the artifacts defined in
[reporting-contract.md](reporting-contract.md). The saved scripts must cover
installation and source verification, download and checksum verification,
capability inspection, profile selection, cohort construction, and analysis.
A Biomni replay alone is not a complete analysis record.
