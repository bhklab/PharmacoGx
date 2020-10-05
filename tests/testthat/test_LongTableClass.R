library(PharmacoGx)
library(parallel)
library(data.table)

## FIXME:: Refactor into multiple test files?

context("Checking LongTable Class Methods.")

# Configure required parameters
filePath <- 'merckLongTable.csv'
rowDataCols <- list(c(cell_line1="cell_line", BatchID="BatchID"))
colDataCols <- list(c(drug1='drugA_name', drug2='drugB_name',
    drug1dose='drugA Conc (µM)', drug2dose='drugB Conc (µM)'))
assayCols <- list(viability=paste0('viability', seq_len(4)),
                  viability_summary=c('mu/muMax', 'X/X0'))
# ---- 1. buildLongTable
test_that('Can build LongTable from single table', {

    context('...from a file path')
    longTable <- buildLongTable(from=filePath, rowDataCols, colDataCols, assayCols)
    expect_s4_class(longTable,'LongTable')
    expect_equal_to_reference(longTable, 'merckLongTable.rds')

    merckDT <- fread(filePath)
    context('...from a data.table')
    longTable1 <- buildLongTable(from=merckDT, rowDataCols, colDataCols, assayCols)
    expect_s4_class(longTable1,'LongTable')
    expect_equal_to_reference(longTable1, 'merckLongTable.rds')

    context('...from a data.frame')
    longTable2 <- buildLongTable(from=setDF(merckDT), rowDataCols, colDataCols, assayCols)
    expect_s4_class(
        longTable2,
        'LongTable')
    expect_equal_to_reference(longTable2, 'merckLongTable.rds')
})

# already know this works from previous test
longTable <- buildLongTable(from=filePath, rowDataCols, colDataCols, assayCols)

# update rowDataCols and colDataCols based on renamed columns in LongTable
rowDataColsOld <- rowDataCols
rowDataCols <- lapply(rowDataCols, names)
colDataColsOlds <- colDataCols
colDataCols <- lapply(colDataCols, names)

## TODO:: Read in raw data for this step so we don't require working accessors
test_that('Can build LongTable from list of tables', {

    assayList <- assays(longTable, withDimnames=TRUE, metadata=TRUE)

    context('...from list of data.table')
    longTable1 <- buildLongTable(from=assayList, rowDataCols, colDataCols, assayCols)
    expect_s4_class(longTable1, 'LongTable')

})

test_that('Function builLongTable errors', {
    context('...assayCols argument missing')
    expect_error(buildLongTable(from=filePath, rowDataCols, colDataCols))

    context('...rowDataCols argument missing')
    expect_error(buildLongTable(from=filePath, colDataCols=colDataCols,
        assayCols=assayCols))

    context('...colDataCols argument missing')
    expect_error(buildLongTable(from=filePath, rowDataCols=rowDataCols,
        assayCols=assayCols))

    context('...assayCols not in data')
    expect_error(buildLongTable(from=filePath, rowDataCols, colDataCols, ))
})

# ---- 2. Acessors methods

test_that('LongTable accessor methods are working', {


})

# ----- 3. Setter Methods

test_that('LongTable setter methods are working', {


})

