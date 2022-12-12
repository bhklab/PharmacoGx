# Package Release News

# 3.3.2
- Debugging vignette issues on the Bioconductor build system

# 3.3.1
- Added new vignette documenting support for drug combination modelling new
drug combination features added in PharmacoGx >=3.0

# 3.1.4
- Modified downloadPSet function to automatically update the PharmacoSet class structure and resave the updated object after download
- This work around is necessary until we can rerun our data engineering pipelines to regenerate all of our PharmacoSet using the 3.1.0 package updates
- Added a number of additional methods for computing drug synergy metrics

# 3.1.1

## 3.1.0
- Update to slot names "cell" -> "sample" and "drug" -> "treatment"
- Update standardized identifier column names to match the above slot nomenclature: "cellid" -> "sampleid", "drugid" -> "treatmentid"

# 2.5

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