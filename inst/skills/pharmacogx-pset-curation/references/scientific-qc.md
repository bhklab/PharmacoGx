# Scientific QC and acceptance evidence

Use the depth of CCLE full-PSet QC and the Nikta/Nasim curation examples as a
standard for exploration and representative analyses when those examples are
available. They are guidance, not API specifications or required folder layouts.
Do not copy their biological answers, thresholds or benchmark results into the
generic workflow. Older examples may bypass accessors or use outdated metric
definitions; verify current package behavior.

## One integrated report

Describe source releases and object contents, then explore actual sample,
treatment, tissue, molecular and response coverage. Show distributions and
missingness, duplicate/mismatched identifiers, excluded rows, replicate support
and fit-failure/censoring rates. Include a few meaningful accessor/analysis
examples chosen before examining significance: response ranking, representative
dose curves, molecular-response matching or combination summaries as supported.
Include ordinary, sparse and problematic records, not only successful fits.

Preserve a compact machine-readable summary behind each figure. Use saved
curated outputs and precomputed validation tables; report rendering should not
download data, refit the entire screen or silently rerun curation. Tests may
check loading, dimensions and keys, but a separate smoke report is not QC.

## Published versus recomputed measurements

Where comparable published data exists:

1. Define the common estimand and match by full source experiment keys,
   including arms, sample, doses, exposure, plate and replicate as relevant.
   Report original counts, matched pairs, exclusions and missingness. Never
   compare by row order or silently average away unmatched dimensions.
2. Establish units, log transforms, normalization, integration range, response
   direction, fit model and threshold definition. Preserve original and
   converted values. If these cannot be reconciled, explain why a numerical
   agreement claim is unavailable.
3. Report correlation **and agreement**: Pearson/Spearman with paired n;
   identity-line scatter; signed bias and absolute error/RMSE in interpretable
   units; and difference-versus-mean plots where appropriate. Correlation alone
   can hide a scale or offset error. Define tolerances from measurement precision
   or scientific use rather than selecting them after seeing the results.
4. Stratify differences by sample/treatment, dose coverage, replicate count and
   fit quality. Inspect outliers and source records; preserve unexplained
   discrepancies instead of tuning methods to maximize agreement.
5. Compare model parameters separately from predictions. Different Hill/slope/
   asymptote parameterizations can yield similar predicted viability over the
   observed range, while similar parameter values can mask unit mistakes.
   Compare fitted predictions on the same dose grid within shared support;
   report prediction error separately from parameter/endpoint error. Mark
   extrapolation and non-identifiable thresholds explicitly.

## Endpoint-specific interpretation

- AUC/AAC require a defined range, axis and normalization. Complementarity is
  conditional on those being identical; do not infer meaning from a column name.
- Check the pinned package's IC50/AC50 definitions: relative half-maximal
  response and absolute viability thresholds are distinct. Source terminology
  may differ. Record out-of-range or unreachable thresholds rather than
  converting them into precise finite potency estimates.
- Preserve dose and viability scales and non-monotonic behavior. Single-dose
  data support observed-response summaries, not an identified dose-response fit.
- For combinations, distinguish observed response, reference expectation,
  difference score and aggregate score. Compare published and recomputed values
  at the same dose-pair and aggregation level; retain monotherapy controls and
  compare multiple relevant reference definitions without selecting a result
  because it is large or significant.
- Molecular QC includes feature/sample mappings, replicate aggregation,
  transformation and assay-specific missingness. Cross-PSet overlap is not
  independent biological replication. A negative predeclared example can pass.

Complete when the artifact reloads, keys/values reconcile, representative
analyses execute with interpretable results, and limitations are documented.
Do not equate `validObject()` success with scientific correctness.
