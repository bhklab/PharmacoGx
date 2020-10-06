# Package Release News

## v2.1.11
- Corrected NEWS.md format as requested by Bioconductor
- Major changes are scheduled for the upcoming Bioconductor release on October 28th, 2020!

## v2.0.0
- PharmacoGx now depends on CoreGx, a package designed to abstract core 
functionality from PharmacoGx for use in other Gx suite packages
- `PharmacoSet` class has been modified to store molecular profile data in the 
`SummarizedExperiment` class instead of the the `ExpressionSet` class
- Argument `pSet` in most PharmacoSet accessor methods now converted to `object` 
instead; this will break code using names parameters for these accessor methods