library(PharmacoGx)
library(qs)

pSetFiles <- list.files('PSets', full.names=TRUE)
pSetNames <- gsub('.*/', '', pSetFiles)
pSetNames <- gsub('.rds', '.qs', pSetFiles)

for (file in pSetFiles) {
    pSet <- readRDS(file)
    newSensitivity <-


}