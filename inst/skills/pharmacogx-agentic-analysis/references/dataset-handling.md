# Dataset-aware PharmacoSet handling

Adapted from Nikta Feizi's updated data skill shared 2026-11-03. These are
release-specific expectations to verify against the selected PSet, not proof
that any particular local object has usable measurements. Preserve original
objects and source values. Use the installed public accessors after checking
their signatures and the object's legacy or TRE/MAE representation.

## Inspect before selecting an analysis

Record the exact PSet release/checksum, biological model, sample and treatment
keys, available response assays/metrics and molecular features. Validate
biological-sample versus assay-sample mappings and complete-case coverage.
Inspect `sampleInfo()`, `treatmentInfo()`, `mDataNames()`, `featureInfo()`,
`phenoInfo()` and `treatmentResponse()` where supported. A failing legacy
`sensitivityMeasures()` accessor does not establish absent response data;
inspect the actual response tables and report the accessor limitation.

Declare units, scale, direction and replicate handling before calculating.
Keep published and recomputed endpoints distinct. In particular, a name
containing `log2` does not establish a log scale; verify source transformations
and distinguish absolute copy number from ploidy-relative measurements.

## Biological model and tissue context

Nikta's canonical-release classifications are:

| Model type | PSet releases |
| --- | --- |
| Established cancer cell lines | CCLE_2015, CTRPv2_2015, GDSC_2020 v1/v2 8.2, PRISM_2020, gCSI_2019, NCI60_2021, FIMM_2016, GRAY_2017, UHNBreast_2019, NCISarcoma_2015, TCL38 |
| Primary AML patient specimens screened ex vivo | BeatAML_2018, Tavor_2020 |
| Xenograft-derived breast tumor cultures (PDTCs) | PDTX_2019 |
| Patient-derived GBM cultures from HGCC | GBM_scr2, GBM_scr3 |

The GBM cultures were established from patient tumors; distinguish them from
freshly isolated patient specimens and conventional established cell-line
panels. Nikta cites Johansson et al., Cell Reports (2020),
https://doi.org/10.1016/j.celrep.2020.107897. Verify biological context for
additional releases rather than extending this table by name alone.

Use consistent, distinguishable colors for model types in cross-dataset plots,
with a legend; use labels, shapes or facets for dataset identity. Report model
type alongside each dataset. Do not label all samples as cell lines or pool
model types without an explicit analysis model and scientific justification.

Consider these disease-specific datasets when they are within the allowed
data scope:

| Context | Candidate PSets |
| --- | --- |
| Breast | GRAY_2017, UHNBreast_2019, PDTX_2019 |
| AML | BeatAML_2018, Tavor_2020 |
| GBM / relevant brain or CNS analyses | GBM_scr2, GBM_scr3 |
| T-cell leukemia or lymphoma | TCL38 |
| Sarcoma | NCISarcoma_2015 |

Also inspect eligible samples from multi-tissue panels. Consider equivalent
tissue terminology and documented single-disease dataset provenance when
annotations are missing. Do not override conflicting sample evidence, treat
pancreas lineage as a PDAC diagnosis, or include TCL38 automatically in AML.
Verify requested drug, assay, feature and final sample coverage for every
candidate. Respect exact cohort definitions in frozen study protocols.

## Response and combination assays

Inspect `treatmentResponse()` and distinguish raw observations, dose-level
synergy, fitted surfaces and context-level `profiles` summaries. Record the
exact assay, score column, model variant and sign. Use dose-level scores for
dose-specific questions. For cell-level associations, derive the prescribed
summary per biological sample and drug pair, accounting for technical contexts
and dose-grid coverage; doses are not independent biological observations.
Use aggregated scores when appropriate to the question or when they are the
only available measurements, with their provenance. Identical metric names
do not guarantee identical scale or computation across PSets.

### GDSC-anchor

Recognize the 7 by 2 anchor-library design. An absent `aac_recomputed_2`
column does not prove the second component has no monotherapy AAC. Inspect
each drug's records and join its corresponding monotherapy measurement to the
same biological sample and relevant experimental context. Preserve both drug
orientations and all tested anchor doses until the analysis defines aggregation.
Use verified Jaaks MSI annotations matched by COSMIC ID when that analysis
requires them; keep missing joins visible. Do not substitute a different
potency-shift definition for a frozen study's estimand.

### GDSC-matrix

Inspect profile column types before association functions. Select compatible
columns into an analysis copy rather than deleting annotations from the source.
Construct dose axes from actual concentrations and preserve drug orientation.

### NCI60 and NCI-ALMANAC

Inspect RNA-seq assay annotations and feature mappings. Correct incompatible
metadata only after confirming the underlying measurement technology, in a
documented copy. Keep source response definitions and modeling transforms
explicit; do not exchange PERCENTGROWTH and PERCENTGROWTHNOTZ or reuse fits
across them. Confirm the requested Bliss/HSA/ZIP variant and dose context.

### BeatAML

Inspect FLT3-ITD sample annotations. Add a derived mutation assay only when
necessary for the requested analysis and supported by those annotations;
preserve missing calls and record derivation. Do not equate missing mutation
information or absence of a reported variant with verified wild type.

### PRISM

Inspect native molecular profiles first. If required measurements are absent,
and external pairing is permitted, follow [PRISM and CCLE](prism-ccle.md).
Record CCLE_2015 explicitly for Nikta's temporary integration route. A newer
CCLE release requires its own validated crosswalk and measurement checks.
Do not treat the original and augmented PRISM objects as independent studies.
If the benchmark asks what PRISM contains natively, retain that native answer;
external augmentation does not change native coverage.

## Eligibility and reporting

Successful matching is not statistical eligibility. Check finite paired
observations, marker/response variability, independent sample counts and group
sizes under the prescribed method. Freeze cohort, endpoint, model, inference,
covariates and multiplicity family before inspecting effects. A frozen
`drugSensitivitySig` ANOVA or tissue-adjusted model takes precedence over the
generic rank-correlation default. Report different valid methods separately
from exact reference concordance, without changing a method to chase agreement.

Report dataset/release/model type, drug/pair, molecular source and feature,
response field, cohort, complete-case and group counts, effect and uncertainty,
raw/adjusted p-values, exclusions and executable provenance. Keep unavailable,
load-failed and unevaluable states distinct from a null biological result.
When screens share CCLE molecular measurements or biological samples, describe
screen reproducibility separately from independent molecular validation.
