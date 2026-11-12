# Analysis decisions and package functions

Read the selected workflow, then inspect the actual PSet and pinned package
help. Examples are illustrative development material; do not tune on held-out
prompts, human scripts or reference results.

## Response summaries and dose curves

Use `sensitivityMeasures()` or the actual TRE assay schema to discover the
endpoint. `summarizeSensitivityProfiles()` can aggregate replicates, but fix
the aggregation and verify output dimnames/orientation. Rank finite supported
records in the correct direction; compare metrics on the same eligible sample
set and report censored values/ties rather than implying interchangeable
estimands. Tiny top-k lists may exhaust the available records.

Before fitting, align dose/viability pairs within sample, treatment, exposure,
batch and replicate. Confirm dose units and `conc_as_log`, viability scale and
`viability_as_pct`, and distinct dose coverage. Inspect raw points and retain
replicate uncertainty. Do not infer a curve from a single dose.

Use `logLogisticRegression()` with a declared `fit_type` and `curve_direction`
where supported. `drugDoseResponseCurve()` can visualize actual/fitted data;
show observations alongside predictions and distinguish interpolation from
extrapolation. Prefer the simplest scientifically adequate fit, not whichever
model maximizes a downstream association.

In the inspected package, `computeAUC()` is normalized viability area and
`computeAAC()` response area, with output scale determined by
`viability_as_pct`. Specify `area.type = "Actual"` or `"Fitted"`, integration
support and truncation. `computeIC50(reference = "relative")` targets relative
half-maximal effect; `computeAC50()` targets absolute 50% viability. Source
publications may use IC50 differently. A fitted threshold outside the measured
range is not well-supported potency. Compare source/recomputed correlation
and agreement only after reconciling units, definitions and support; compare
curve predictions separately from fitted parameters.

## Biomarkers and cross-PSet replication

Distinguish continuous expression/CNV from binary or categorical mutation,
fusion and other events. Do not use mutation as a substitute for fusion, or
gene-level aggregation as proof that a particular variant is present.
Audit sample mapping, tissue/disease context, finite pairs and group support.
Use the method guidance in [direct-r.md](direct-r.md); metadata associations
are not molecular biomarkers. Define candidate features and multiplicity
before inspecting results and label broad screens as exploratory.

Use `intersectPSet()` or explicit identifier mappings only after verifying
identity and endpoint compatibility. Preserve each dataset's original sample,
treatment and feature IDs. Analyze separately and compare signed effects,
uncertainty and sample overlap. Repeated cell lines in another screen are
technical/contextual replication, not independent biological validation.
Do not pool raw values or select discovery/validation roles after seeing results.

## Combinations and mono-versus-combo

Inspect `treatmentResponse()` to distinguish raw response, mono profiles,
combo profiles and published/recomputed synergy assays. Retain arm order or
an explicit symmetry rule and match controls by dose, time, plate and replicate
when possible. Dose-pair scores and pair-level summaries are different objects.

`computeBliss()` and `computeHSA()` return reference viability, not an observed
synergy score. For fractional viability, Bliss multiplies matched monotherapy
viabilities and HSA takes their elementwise minimum. A declared
reference-minus-observed score is positive for extra killing. Never multiply
percentages as though they were fractions. Loewe/ZIP need their documented
inputs, fitted monotherapies and unit conventions; inspect `computeLoewe()`,
`computeZIP()`/`computeZIPdelta()` help rather than inventing missing curves.

Rank pairs using prespecified sample support, robust summaries and consistency;
show uncertainty and alternatives instead of the largest single score. A
selected heatmap is one sample's view. For combination-specific biomarkers,
use matched samples and compare the combo phenotype with both monotherapy arms.
Differences in significance alone are not evidence that effects differ; define
the interaction/contrast needed for the claim and account for multiple testing.

## Paper-study reproduction

Build a discrepancy ledger for source release, cohort, feature definition,
response endpoint, preprocessing, model/covariates, inference and multiplicity.
Freeze the faithful primary reproduction; label any modernization or alternate
method separately. Compare numerical results only when estimands and matching
agree; distinguish exact reproduction, method adaptation and biological
replication. A negative or unavailable result is informative, not permission
to search for a replacement significant example.

For prediction tasks, split by the biological unit before preprocessing and
feature selection; preserve an independent test set and avoid sample/drug
overlap leakage. Association discovery is not out-of-sample predictive accuracy.
Keep evaluator-only reference scripts/results outside participant resources.
