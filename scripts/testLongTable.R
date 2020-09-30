library(PharmacoGx)
library(data.table)

rowDataCols <- list(c(cell_line1="cell_line", BatchID="BatchID"))
colDataCols <- list(c(drug1='drugA_name', drug2='drugB_name', drug1dose='drugA Conc (µM)', drug2dose='drugB Conc (µM)'))
filePath <- 'data/drug_combo_merck.csv'
assayCols <- list(viability=paste0('viability', seq_len(4)),
                  viability_summary=c('mu/muMax', 'X/X0'))

a <- Sys.time()
longTable <- buildLongTable(filePath, rowDataCols, colDataCols, assayCols)
b <- Sys.time()
b - a

a <- Sys.time()
reindex(longTable)  # currently a tiny bit slower than building a long table
b <- Sys.time()
b - a

longTable['VCAP', ]

longTable[.(cell_line == 'VCAP' & BatchID == 1),
          .(drugA_name == "Lapatinib" & !(drugB_name == 'Dasatinib'))]

longTable[['viability']]

longTable$viability[, colMeans(.SD, na.rm=TRUE)]
