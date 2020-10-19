library(PharmacoGx)
library(qs)
library(data.table)

setDTthreads(15)

pSetFiles <- list.files('PSets', pattern='.*rds', full.names=TRUE)
pSetNames <- gsub('.*/', '', pSetFiles)
pSetNames <- gsub('.rds', '.qs', pSetNames)
i <- 1

#for (i in seq_along(pSetFiles)) {
#    print(pSetNames[i])
#    pSet <- readRDS(pSetFiles[i])
#    .pSet <- pSet
#    newSensitivity <- .sensitivitySlotToLongTable(pSet)
#    pSet@sensitivity <- newSensitivity
#    qsave(pSet, file.path('PSets', pSetNames[i]), nthreads=15)
#}

pSet <- readRDS(pSetFiles[1])
.pSet <- pSet

LT <- sensitivitySlotToLongTable(pSet)
pSet@sensitivity <- LT

OLT <- pSet@sensitivity

oldRaw <- sensitivityRaw(.pSet)
sensitivityRaw(pSet) <- oldRaw
raw <- sensitivityRaw(pSet)

all.equal(raw[rownames(oldRaw),,], oldRaw)








class(sensitivitySlot(.pSet))
oldProfSum <- summarizeSensitivityProfiles(.pSet, 'aac_recomputed', summary.stat='mean')
class(sensitivitySlot(pSet))
profSum <- summarizeSensitivityProfiles(pSet, 'aac_recomputed', summary.stat='mean')

all.equal(oldProfSum, profSum[rownames(oldProfSum), colnames(oldProfSum)])



prof <- sensitivityProfiles(pSet)
sensitivityProfiles(pSet) <- prof
SP <- sensitivityProfiles(pSet)
all.equal(SP[rownames(prof), colnames(prof)], prof)

LT <- pSet@sensitivity

for (i in 1:100) {
    n <- sample(1:nrow(prof), 1)
    print(prof[n, ])
    print(SP[rownames(prof)[n], ])
    Sys.sleep(1)
}



LTold <- sensitivitySlot(.pSet)
LT <- sensitivitySlot(pSet)

pSet <- .pSet
info <- sensitivityInfo(pSet)

sensitivityInfo(pSet) <- info
SI <- sensitivityInfo(pSet)

all.equal(SI[rownames(info), colnames(info)], info)
