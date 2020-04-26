# Package Release News

## v2.0.0
- PharmacoGx now depends on CoreGx, a package designed to abstract core 
functionality from PharmacoGx for use in other Gx suite packages
- `PharmacoSet` class has been modified to store molecular profile data in the 
`SummarizedExperiment` class instead of the the `ExpressionSet` class
- Argument `pSet` in most PharmacoSet accessor methods now converted to `object` 
instead; this will break code using names parameters for these accessor methods

## Future Releases
- Convert the `@molecularProfiles` slot to be a `MultiAssayExperiment` object 
instead of a `list`
- Reimplement `PharmacoSet` set operations (intersect, subset, etc.) using 
generics defined in the `BiocGenerics` package
- Redesign the `@sensitivity` and `@perturbation` slots to use either existing
Bioconductor classes or custom implementations of Bioconductor classes
- Add support for drug combinations to the package
- Add support for new types of molecular data such as metabolomics, proteomics,
and single cell sequencing
- Continue to abstract functionality into CoreGx