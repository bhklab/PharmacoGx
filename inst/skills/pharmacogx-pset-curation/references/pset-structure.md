# PSet structure and version-aware R access

The API notes below were checked against PharmacoGx `feat/mcp` source
`84e8099fa508805f2f62d1e618458c9b78364ca5` and CoreGx 2.16.0 on 2026-10-27.
Bounded constructor, accessor and source-value round-trip checks also passed
using installed PR178 PharmacoGx `28c58707eece80086bee75d8ef0e8afae96389ca`
and CoreGx 2.7.3 `bf6d814abd3dc0fc98f9f419f430cc6473046cb2`, with loaded
namespace paths and `RemoteSha` verified in an isolated R 4.5.3 runtime.
Recheck signatures for the environment actually used. The package's
`R/PharmacoSet-class.R`, `R/methods-summarizeSensitivityProfiles.R` and CoreGx's
`TreatmentResponseExperiment` help are the implementation references.

## Object boundaries

- `PharmacoSet()` supports the legacy sensitivityInfo/sensitivityRaw/
  sensitivityProfiles constructor path. `PharmacoSet2()` accepts explicit
  `sample`, `treatment`, `molecularProfiles`, `treatmentResponse`, `curation`
  and `datasetType`. Both return a PharmacoSet; do not infer storage from class
  name alone or invent a `PharmacoSet3()` constructor.
- The modern path uses `MultiAssayExperiment` molecular profiles containing
  `SummarizedExperiment` assays and a `CoreGx::TreatmentResponseExperiment`
  (TRE), currently built on LongTable. A TRE-only object is not a full PSet.
- Molecular assay rows are features and columns are assay samples. Align
  `rowData`/`colData` exactly; use the MAE `sampleMap` to relate assay columns
  to primary samples. Check legacy `cellid`/`sampleid` annotations explicitly.
  Resolve conflicts rather than silently taking the first ID field.
- The tested CoreGx 2.16.0 molecular-profile validator requires `colData$batchid`.
  Map a supplied batch identifier explicitly (for example `batch_id` to
  `batchid`); do not invent a biological batch to satisfy validation. If the
  source lacks it, document the missing metadata and investigate the pinned
  constructor's supported representation.
- `metadata(se)$annotation` is interpreted by some PharmacoGx molecular
  methods. Use the vocabulary supported by the pinned package, consistent with
  the real measurement; a profile name is not evidence of usable measurements.
- TRE `rowIDs`, `colIDs` and per-assay `assayIDs` declare keys. The row/column
  orientation is configured, not a universal rectangular drug-by-cell matrix.
  Put entity attributes in metadata tables; assay tables hold keys and values.
  An assay key must include the dimensions needed to distinguish each record,
  including dose, time, plate or replicate when these vary.

Both tested CoreGx builds need particular care: a bounded constructor test with dose and
replicate dimensions only in `assayIDs` collapsed 12 source rows to two, while
`validObject()` still passed. In both environments, putting dose/time in
`rowIDs` and replicate in `colIDs`, with corresponding metadata rows, preserved
all measurements. Choose the partition appropriate to the data and pinned
implementation; keep primary sample/treatment IDs separate from measurement
keys. Verify every source key and value after construction and save/reload.
Do not assume a newer source checkout describes the installed constructor.

## Inspection recipe

Use bounded previews; do not print a full molecular matrix.

```r
library(PharmacoGx)
stopifnot(methods::is(pset, "PharmacoSet"))
methods::validObject(pset)
samples <- sampleInfo(pset)
treatments <- treatmentInfo(pset)
profiles <- mDataNames(pset)
response <- treatmentResponse(pset)
list(samples = dim(samples), treatments = dim(treatments),
     molecular_profiles = profiles, response_class = class(response))
```

For a TRE, inspect available assays and their key metadata before extraction;
use CoreGx/SummarizedExperiment accessor methods supported by that version.
For legacy storage, inspect `sensitivityInfo()`, `sensitivityMeasures()` and
`sensitivityProfiles()`. `summarizeSensitivityProfiles()` supports LongTable
with explicit `profiles_assay`, `treatment_col`, `sample_col` and
`sensitivity.measure`; verify returned dimnames rather than assuming matrix
orientation. Combination rows must not silently become monotherapy summaries.

## Constructor shape

This is the parameter contract, not a complete dataset-specific pipeline:

```r
tre <- CoreGx::TreatmentResponseExperiment(
  rowData = treatment_keys, rowIDs = treatment_key_columns,
  colData = sample_keys, colIDs = sample_key_columns,
  assays = response_tables, assayIDs = response_key_columns,
  metadata = response_metadata
)
pset <- PharmacoGx::PharmacoSet2(
  name = dataset_name, sample = sample_metadata,
  treatment = treatment_metadata, molecularProfiles = molecular_mae,
  treatmentResponse = tre,
  curation = list(sample = sample_mapping, treatment = treatment_mapping,
                  tissue = tissue_mapping),
  datasetType = "sensitivity"
)
```

Check all treatment arms against treatment metadata and all primary samples
against sample metadata; constructor success may not enforce every scientific
relationship. For each assay compare source and curated counts, key sets,
missingness and selected values. Test one-to-many joins before expanding data;
an unintended Cartesian join can produce a structurally valid wrong object.

Use `updateObject()` only as a documented migration into a new output; retain
the original bytes/checksum. Current `downloadPSet()` loads, updates and rewrites
its file. To preserve source-byte provenance, acquire and checksum an immutable
source copy separately before constructing a migrated artifact.
