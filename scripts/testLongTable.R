library(PharmacoGx)
library(data.table)
library(microbenchmark)

setDTthreads(2)

rowDataCols1 <- list(c(cell_line1="cell_line", BatchID="BatchID"))
colDataCols1 <- list(c(drug1='drugA_name', drug2='drugB_name',
    drug1dose='drugA Conc (µM)', drug2dose='drugB Conc (µM)'))
filePath <- 'data/drug_combo_merck.csv'
assayCols <- list(#dose=c('drugA Conc (µM)', 'drugB Conc (µM)'),
                  viability=paste0('viability', seq_len(4)),
                  viability_summary=c('mu/muMax', 'X/X0'))

microbenchmark({
    longTable <- buildLongTable(filePath, rowDataCols1, colDataCols1, assayCols)
}, times=1L)

# rowdata <- rowData(longTable)
# rowdata[, testCol := rnorm(nrow(longTable))]
# rowdata
#
# rowData(longTable) <- rowdata
# rowData(longTable)
#
# coldata <- colData(longTable)
# coldata[, testCols := rnorm(ncol(longTable))]
# coldata

# colData(longTable) <- coldata
# colData(longTable)

## TODO:: This needs a lot of parameters now?
from <- assays(longTable, withDimnames=TRUE, metadata=TRUE, key=FALSE)
from$new_viab <- from$viability

assayCols$new_viab <- assayCols(longTable, 'viability')


## TODO:: Assay is slower than assays for assignment
#a <- Sys.time()
## FIXME:: This is adding prefixes to assay column names?
#assay(longTable, 'new_viability') <- viab
#b <- Sys.time()
#b - a

rowDataCols <- lapply(rowDataCols1, names)
colDataCols <- lapply(colDataCols1, names)
microbenchmark({
    longTable2 <- buildLongTable(from, rowDataCols, colDataCols, assayCols)
}, times=10L)


a <- Sys.time()
assays(longTable) <- from
b <- Sys.time()
b - a

a <- Sys.time()
reindex(longTable)  # currently a tiny bit slower than building a long table
b <- Sys.time()
b - a

longTable['VCAP', ]

longTable[.(cell_line1 == 'VCAP' & BatchID == 1),
          .(drug1 == "Lapatinib" & !(drug2 == 'Dasatinib'))]

longTable[['viability']]

longTable[, c('5-FU:*', 'Temozolamide:MK-8776')]

x <- longTable$viability[, mean := rowMeans(.SD, na.rm=TRUE)]

# Compare runtimes of subset with and without reindexing
a <- Sys.time()
subsetTo(longTable, .(cell_line1 == 'VCAP' & BatchID != 2), .(drug1 != 'Lapatinib', reindex=FALSE))
b <- Sys.time()
t1 <- b - a

a <- Sys.time()
subsetTo(longTable, .(cell_line1 == 'VCAP' & BatchID != 2), .(drug1 != 'Lapatinib', reindex=TRUE))
b <- Sys.time()
t2 <- b - a


as.numeric(t2)/as.numeric(t1)