library(PharmacoGx)
library(data.table)

rowDataCols <- list(c("cell_line"), c("BatchID"))
colDataCols <- list(c("drugA_name", "drugB_name"))
filePath <- 'data/drug_combo_merck.csv'
assayCols <- list(dose=c('drugA Conc (µM)', 'drugB Conc (µM)'),
                  viability=paste0('viability', seq_len(4)),
                  viability_summary=c('mu/muMax', 'X/X0'))

longTable <- buildLongTableFromCSV(filePath, rowDataCols, colDataCols, assayCols)


longTable['VCAP', ]

longTable[['viability']]

longTable$viability[, colMeans(.SD, na.rm=TRUE)]
