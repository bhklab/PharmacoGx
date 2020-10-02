library(PharmacoGx)
library(data.table)

rowDataCols1 <- list(c(cell_line1="cell_line", BatchID="BatchID"))
colDataCols1 <- list(c(drug1='drugA_name', drug2='drugB_name',
    drug1dose='drugA Conc (µM)', drug2dose='drugB Conc (µM)'))
filePath <- 'data/drug_combo_merck.csv'
assayCols <- list(viability=paste0('viability', seq_len(4)),
                  viability_summary=c('mu/muMax', 'X/X0'))

a <- Sys.time()
longTable <- buildLongTable(filePath, rowDataCols1, colDataCols1, assayCols)
b <- Sys.time()
b - a

viab <- assay(longTable, 'viability', withDimnames=TRUE, metadata=TRUE)

## TODO:: Assay is slower than assays for assignment
a <- Sys.time()
assay(longTable, 'new_viability') <- viab
b <- Sys.time()
b - a

from <- assays(longTable, withDimnames=TRUE, metadata=TRUE)

assayCols$new_viability <- assayCols$viability

rowDataCols <- lapply(rowDataCols1, names)
colDataCols <- lapply(colDataCols1, names)
longTable2 <- buildLongTable(from, rowDataCols, colDataCols, assayCols)

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
subset(longTable, .(cell_line1 == 'VCAP' & BatchID != 2), .(drug1 != 'Lapatinib', reindex=FALSE))
b <- Sys.time()
t1 <- b - a

a <- Sys.time()
subset(longTable, .(cell_line1 == 'VCAP' & BatchID != 2), .(drug1 != 'Lapatinib', reindex=TRUE))
b <- Sys.time()
t2 <- b - a


as.numeric(t2)/as.numeric(t1)