library(PharmacoGx)
library(parallel)
library(data.table)

## FIXME:: Refactor into multiple test files?

context('Checking LongTable Class Methods.')

# Configure required parameters
filePath <- 'merckLongTable.csv'
rowDataCols <- list(c(cell_line1='cell_line', BatchID='BatchID'))
colDataCols <- list(c(drug1='drugA_name', drug2='drugB_name',
    drug1dose='drugA Conc (µM)', drug2dose='drugB Conc (µM)'))
assayCols <- list(viability=paste0('viability', seq_len(4)),
                  viability_summary=c('mu/muMax', 'X/X0'))
# ---- 1. buildLongTable
context('Testing buildLongTable function')

context('   Testing buildLongTable from a single table')
test_that('Can build LongTable from single table', {

    test_that('    from a file path', {
        longTable <- buildLongTable(from=filePath, rowDataCols, colDataCols, assayCols)
        expect_s4_class(longTable,'LongTable')
        expect_equal_to_reference(longTable, 'merckLongTable.rds')
    })

    merckDT <- fread(filePath)
    test_that('   from a data.table', {
        longTable1 <- buildLongTable(from=merckDT, rowDataCols, colDataCols, assayCols)
        expect_s4_class(longTable1,'LongTable')
        expect_equal_to_reference(longTable1, 'merckLongTable.rds')
    })


    test_that('   from a data.frame', {
        longTable2 <- buildLongTable(from=setDF(merckDT), rowDataCols, colDataCols, assayCols)
        expect_s4_class(longTable2, 'LongTable')
        expect_equal_to_reference(longTable2, 'merckLongTable.rds')
    })
})

# already know this works from previous test
longTable <- buildLongTable(from=filePath, rowDataCols, colDataCols, assayCols)

# update rowDataCols and colDataCols based on renamed columns in LongTable
rowDataColsOld <- rowDataCols
rowDataCols <- lapply(rowDataCols, names)
colDataColsOlds <- colDataCols
colDataCols <- lapply(colDataCols, names)

## TODO:: Read in raw data for this step so we don't require working accessors
context('    Testing buildLongTable from mutliple tables')
test_that('Can build LongTable from list of tables', {

    assayList <- assays(longTable, withDimnames=TRUE, metadata=TRUE)

    test_that('    from list of data.table', {
        longTable1 <- buildLongTable(from=assayList, rowDataCols, colDataCols, assayCols)
        expect_s4_class(longTable1, 'LongTable')
        expect_equal_to_reference(longTable1, 'merckLongTable.rds')
    })

})

context('   buildLongTable errors')
test_that('Function buildLongTable errors', {
    test_that('    assayCols argument missing', {
        expect_error(buildLongTable(from=filePath, rowDataCols, colDataCols))
    })

    test_that('    rowDataCols argument missing', {
        expect_error(buildLongTable(from=filePath, colDataCols=colDataCols,
            assayCols=assayCols))
    })

    test_that('    colDataCols argument missing', {
        expect_error(buildLongTable(from=filePath, rowDataCols=rowDataCols,
            assayCols=assayCols))
    })

    test_that('    assayCols not in data', {
        expect_error(buildLongTable(from=filePath, rowDataCols, colDataCols, ))
    })

    ## TODO:: Add more informative error cases

})

# ---- 2. Acessors methods

# work correctly
context('Testing LongTable accessor methods')

context('   LongTable accessors work correctly')
test_that('LongTable accessor methods are working', {

    test_that('rowData works correctly', {
        # without keys
        rowData <- rowData(longTable)
        expect_equal(colnames(rowData), c('cell_line1', 'BatchID'))
        expect_s3_class(rowData, 'data.table')
        expect_equal_to_reference(rowData, 'merckLongTable.rowData.rds')

        # with keys
        rowDataKey <- rowData(longTable, key=TRUE)
        expect_equal(colnames(rowDataKey), c('cell_line1', 'BatchID', 'rowKey'))
        expect_equal_to_reference(rowDataKey, 'merckLongTable.rowDataWithKey.rds')
    })

    test_that('colData works correctly', {
        # without keys
        colData <- colData(longTable)
        expect_equal(colnames(colData), c("drug1", "drug2", "drug1dose", "drug2dose"))
        expect_s3_class(colData, 'data.table')
        expect_equal_to_reference(colData, 'merckLongTable.colData.rds')

        # with keys
        colDataKey <- colData(longTable, key=TRUE)
        expect_equal(colnames(colDataKey), c("drug1", "drug2", "drug1dose", "drug2dose", 'colKey'))
        expect_equal_to_reference(colDataKey, 'merckLongTable.colDataWithKey.rds')
    })

    test_that('assays works correctly', {
        assays <- assays(longTable)
        expect_equal(names(assays), assayNames(longTable))
        expect_equal_to_reference(assays, 'merckLongTable.assays.rds')

        assaysDimnames <- assays(longTable, withDimnames=TRUE)
        expect_equal_to_reference(assaysDimnames, 'merckLongTable.assaysDimnames.rds')

        assaysDimnamesMetadata <- assays(longTable, withDimnames=TRUE, metadata=TRUE)
        expect_equal_to_reference(assaysDimnamesMetadata, 'merckLongTable.assayDimnamesMetdata.rds')
    })

    test_that('assay works correctly', {
        assay <- assay(longTable, assayNames(longTable)[1])
        expect_equal_to_reference(assay, 'merckLongTable.assay.viability.rds')

        assayDimnames <- assay(longTable, assayNames(longTable)[1], withDimnames=TRUE)
        expect_equal_to_reference(assayDimnames, 'merckLongTable.assayDimnames.viability.rds')
    })

})

# error correctly
context('   LongTable accessor methods error')

# ----- 3. Setter Methods

# work correctly
context('Testing LongTable setter methods')

test_that('LongTable setter methods are working', {


})

# error correctly
context('    LongTable setter methods error')

test_that('LongTable accessor methods error correctly', {

})

# ---- 4. Subset methods