# Package Release News

# 3.0.1
- Fix bug in logLogsticRegression function causing tests to fail due to loss of numerical precision in floating point values

# 3.0.0
- Bioconductor fall 2022 release!
- Includes updated TreatmentResponeExperiment class for storing drug combination experiments.
## 2.5.3
- Added PharmacoSet2 constructor to allow creation of PSets with updated class definition introducted in BioC 3.13
- The sensitivity slot is now required to be a TreatmentResponseExperiment
- The molecularProfiles slot is now required to be a MultiAssayExperiment
- The original constructor and all accessors remain in the package for full backwards compatibility

## v2.5.2
- Fix: remove 'fdr' item from geneDrugSensitivity return vector

## v2.5.1
- Fix: reverted GDSCsmall.rda and CCLEsmall.rda to original data format; they
were accidentally pushed with MultiAssayExperiments in @molecularProfiles

## v2.5.0
- Spring Bioconductor release!

## v2.1.12
- Added experimental support for a new class, the `LongTable`, for storing the
sensitivity data in a `PharmacoSet`.
- Because this is not well tested, we have left not updated the PSets available
via the `downloadPSets` function.
- Instead we have provided a convenient function,
`convertSensitivitySlotToLongTable`, which takes in a `PharmacoSet` object,
converts the data in the `@sensitivty` slot to a `LongTable` and returns an
updated `PharmacoSet` object.
- The `LongTable` class will be used in the future to allow `PharmacoSet`s to
store treatment response experiments with multiple drugs or cell-lines, greatly
expanding the variety of data which can be stored in a `PharmacoSet` object.
- For more details on the `LongTable` class, please see the vignette in the
`CoreGx` package.

## v2.0.0
- PharmacoGx now depends on CoreGx, a package designed to abstract core
functionality from PharmacoGx for use in other Gx suite packages
- `PharmacoSet` class has been modified to store molecular profile data in the
`SummarizedExperiment` class instead of the the `ExpressionSet` class
- Argument `pSet` in most PharmacoSet accessor methods now converted to `object`
instead; this will break code using names parameters for these accessor methods