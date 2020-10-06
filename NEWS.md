# Package Release News

## v2.0.0
- PharmacoGx now depends on CoreGx, a package designed to abstract core 
functionality from PharmacoGx for use in other Gx suite packages
- `PharmacoSet` class has been modified to store molecular profile data in the 
`SummarizedExperiment` class instead of the the `ExpressionSet` class
- Argument `pSet` in most PharmacoSet accessor methods now converted to `object` 
instead; this will break code using names parameters for these accessor methods