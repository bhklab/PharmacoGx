library(PharmacoGx)
library(data.table)

rowDataCols <- list(c("cell_line"), c("BatchID"))
colDataCols <- list(c("drugA_name", "drugB_name"))
filePath <- 'data/drug_combo_merck.csv'
assayCols <- list(dose=c('drugA Conc (µM)', 'drugB Conc (µM)'),
                  viability=paste0('viability', seq_len(4)),
                  viability_summary=c('mu/muMax', 'X/X0'))

a <- Sys.time()
longTable <- buildLongTableFromCSV(filePath, rowDataCols, colDataCols, assayCols)
b <- Sys.time()
b - a

a <- Sys.time()
reindex(longTable)  # currently a tiny bit slower than building a long table
b <- Sys.time()
b - a

longTable['VCAP', ]

longTable[.(cellLine1 == 'VCAP' & BatchID == 1),
          .(drug1 == "Lapatinib" & !(drug2 == 'Dasatinib'))]

longTable[['viability']]

longTable$viability[, colMeans(.SD, na.rm=TRUE)]
