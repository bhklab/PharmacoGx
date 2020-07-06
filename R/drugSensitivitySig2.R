.drugSensitivitySig2PharmacoSet <- function(object, mDataType) {

    .data.table.rn <- function(...) data.table(..., keep.rownames='rownames')

    sensProf <-  .data.table.rn(sensitivityProfiles(object))
    sensInf <-   .data.table.rn(sensitivityInfo(object))

    sensData <- merge(sensInf, sensProf, on="rownames")

    cellInf <- .data.table.rn(cellInfo(object))

    sensDT <- merge(sensData, cellInf, by.x='cellid', by.y='cell_id')

    mProf <- .data.table.rn(molecularProfiles(object, mDataType))
    fInf <- data.table(as(featureInfo(object, mDataType), 'data.frame'))
    pInf <- data.table(as(phenoInfo(object, mDataType), 'data.frame'))

    .melt <- function(...) data.table::melt(..., variable.factor=FALSE)

    longMolProf <- .melt(mProf, measure.vars=colnames(mProf)[2:ncol(mProf)], variable.name='sample_id', value.name="expression")
    longMolProfFinfo <- merge(longMolProf, fInf, by='rownames', allow.cartesian=TRUE)
    longSE <- merge(longMolProfFinfo, pInf, by.x='sample_id', by.y="rownames")

    longDT <- merge(longSE, sensDT[unique(cellid)], by="cellid", allow.cartesian=TRUE)

    DT <- longDT[, .(cellid, sample_id, expression, EnsemblGeneId, batchid, drug.name, auc_recomputed, tissueid)]

    permutationsTest <- function(someInput) { "Do some calculation" }

    drugSensitivity <- DT[, lapply(.SD, perturmutationsTest), by=.(EnsemblGeneId, drug.name, tissueid)]
}

.summarizedExperimentToDT <- function(SE) {
    mProf <- .data.table.rn(molecularProfiles(object, mDataType))
    fInf <- data.table(as(featureInfo(object, mDataType), 'data.frame'))
    pInf <- data.table(as(phenoInfo(object, mDataType), 'data.frame'))

    .melt <- function(...) data.table::melt(..., variable.factor=FALSE)

    longMolProf <- .melt(mProf, measure.vars=colnames(mProf)[2:ncol(mProf)], variable.name='sample_id', value.name="expression")
    longMolProfFinfo <- merge(longMolProf, fInf, by='rownames', allow.cartesian=TRUE)
    longSE <- merge(longMolProfFinfo, pInf, by.x='sample_id', by.y="rownames")

    longDT <- merge(longSE, sensDT[unique(cellid)], by="cellid", allow.cartesian=TRUE)
}