library(PharmacoGx)
library(qs)

pSetFiles <- list.files('PSets', full.names=TRUE)
pSetNames <- gsub('.*/', '', pSetFiles)
pSetNames <- gsub('.rds', '.qs', pSetFiles)

for (i in seq_along(pSetFiles)) {
    print(pSetNames[i])
    pSet <- readRDS(pSetFiles[i])
    newSensitivity <- .sensitivitySlotToLongTable(pSet)
    pSet@sensitivity <- newSensitivity
    qsave(pSet, file.path('NewPSets', pSetNames[i]))
}