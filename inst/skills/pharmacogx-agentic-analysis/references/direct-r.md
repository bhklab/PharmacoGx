# Direct R in Claude, Biomni and other hosts

Use an R-capable execution environment. Loading a Skill does not install R,
packages, datasets or a connector. Begin with a fresh R process and a dedicated
output directory; retain the full script, command, warnings, errors and exit
status. Use `Rscript --vanilla analysis.R` when a shell is available. Start
with one thread and load large PSets sequentially; release objects or processes
between datasets. Do not change a running curation environment.

Some restricted runtimes return `NA` from `parallel::detectCores()`. The tested
`drugSensitivitySig()` implementation then fails even with `nthread = 1`.
Treat this as an execution-environment failure; use an approved runtime with
working CPU detection or report the blocked analysis. Do not silently replace
the statistical method to obtain a result.

## Verify the requested package source

Record `R.version.string`, installed package paths, `packageVersion()` for
PharmacoGx/CoreGx and the dependency lock or `sessionInfo()`. For Git installs,
record DESCRIPTION's `RemoteUsername`, `RemoteRepo`, `RemoteRef`, `RemoteSha`.
Compare the SHA with the exact task-approved revision. A moving branch name,
release version or unverified preinstalled package is insufficient when the task
requires development features. Do not silently substitute a release package.

For an approved source checkout, record its full Git SHA and dirty status,
then explicitly load that checkout with `pkgload::load_all()` in the isolated
process. This is a separate source-verification path; missing `RemoteSha` in a
checkout does not by itself mean failure. For installation, use `pak` with the
approved SHA in a designated environment; do not install or upgrade into an
active shared analysis environment. Restart R after installation.

Never source executable R from an unreviewed remote URL. Inspect documentation
and method signatures from the verified package before adapting older scripts.

On 2026-10-27 the bounded `scripts/direct-r-examples.R` checks passed in an
isolated R 4.5.3 runtime with loaded PharmacoGx PR178
`28c58707eece80086bee75d8ef0e8afae96389ca` and CoreGx
`bf6d814abd3dc0fc98f9f419f430cc6473046cb2` (2.7.3), verified through namespace
paths and installed `RemoteSha`. This covers the bundled response/signature and
synthetic area, curve, potency and combination-reference examples; it does not
establish scientific validity for a new PSet or certify every package API.

## Inspect, then analyze

Use `sampleInfo()`, `treatmentInfo()`, `mDataNames()`, `featureInfo()`,
`molecularProfiles()` and `treatmentResponse()` to inspect the object. Legacy
sensitivity accessors and modern TRE tables differ. Validate assay IDs against
biological sample IDs via annotations/sampleMap; a zero-column extraction can
be a mapping failure rather than absent biology. Do not silently borrow an MCP
helper and call it validated direct-package execution.

Resolve exact treatment/feature IDs, assay units and transformations, disease
annotations and the response endpoint. Enumerate relevant molecular profiles,
their cohort coverage and finite complete pairs. Select using a rule fixed
before association results, not the first profile or best p-value. Save
`profile_selection.csv` and freeze the analysis manifest.

`availablePSets()` is a distribution catalog, not a joint coverage index.
Remote discovery metadata requires local PSet verification. Use existing
verified files and staged approved downloads; preserve original download bytes
and checksums because `downloadPSet()` can migrate and rewrite its file.

## Molecular associations

For `drugSensitivitySig()` explicitly fix `mDataType`, `drugs`, `features`,
`cells`/`tissues`, `sensitivity.measure`, `molecular.summary.stat`,
`sensitivity.summary.stat`, `standardize`, `modeling.method`,
`inference.method`, `returnValues` and `nthread`. Preserve covariates and record
the complete call. Use named result dimensions such as `estimate`, `pvalue`,
`fdr`, `lower`, `upper` or `se` only when actually returned; never index the
internal `.Data` slot by position.

Its `modeling.method = "spearman"` uses ranked inputs in a covariate-adjusted
partial-correlation workflow, not automatically simple pairwise Spearman.
`"anova"`/`"lm"` compare nested models; this is not automatically the same as
every published one-way ANOVA. `inference.method = "analytic"` and
`"resampling"` are distinct; check supported combinations in the current API.
Specify `standardize` rather than letting it change the estimand unnoticed.

Use `stats::cor.test(..., method = "spearman", exact = FALSE)` on explicitly
matched finite pairs when simple asymptotic Spearman is the frozen method.
For Welch comparisons use explicit `stats::t.test(..., var.equal = FALSE)`
with declared independent groups. For a requested Glass's delta, state which
group supplies the standard deviation and the signed numerator; report group
sizes and uncertainty and handle zero SD. Match the paper's method rather than
substituting a convenient MCP association tool.

Analyze each PSet independently; define adjustment across the complete
confirmatory family, including separate calls/datasets where appropriate.
Report uncertainty, attrition, effect direction and weak/negative findings.
Follow [reporting-contract.md](reporting-contract.md).
