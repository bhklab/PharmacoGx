pkgname <- "PharmacoGx"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('PharmacoGx')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("PharmacoSet")
### * PharmacoSet

flush(stderr()); flush(stdout())

### Name: PharmacoSet
### Title: PharmacoSet constructor
### Aliases: PharmacoSet

### ** Examples

 
## For help creating a PharmacoSet object, please see the following vignette:
browseVignettes("PharmacoGx")




cleanEx()
nameEx("amcc")
### * amcc

flush(stderr()); flush(stdout())

### Name: amcc
### Title: Calculate an Adaptive Matthews Correlation Coefficient
### Aliases: amcc

### ** Examples

x <- c(1,2,3,4,5,6,7)
y <- c(1,3,5,4,2,7,6)
amcc(x,y, min.cat=2)




cleanEx()
nameEx("availablePSets")
### * availablePSets

flush(stderr()); flush(stdout())

### Name: availablePSets
### Title: Return a table of PharmacoSets available for download
### Aliases: availablePSets

### ** Examples

if (interactive()){
availablePSets()
}




cleanEx()
nameEx("cellInfo-set")
### * cellInfo-set

flush(stderr()); flush(stdout())

### Name: cellInfo<-
### Title: cellInfo<- Generic
### Aliases: cellInfo<-

### ** Examples

data(CCLEsmall)
cellInfo(CCLEsmall) <- cellInfo(CCLEsmall)




cleanEx()
nameEx("cellInfo")
### * cellInfo

flush(stderr()); flush(stdout())

### Name: cellInfo
### Title: cellInfo Generic
### Aliases: cellInfo

### ** Examples

data(CCLEsmall)
cellInfo(CCLEsmall)




cleanEx()
nameEx("cellNames-set")
### * cellNames-set

flush(stderr()); flush(stdout())

### Name: cellNames<-
### Title: cellNames<- Generic
### Aliases: cellNames<-

### ** Examples

data(CCLEsmall)
cellNames(CCLEsmall) <- cellNames(CCLEsmall)




cleanEx()
nameEx("cellNames")
### * cellNames

flush(stderr()); flush(stdout())

### Name: cellNames
### Title: cellNames Generic
### Aliases: cellNames

### ** Examples

data(CCLEsmall)
cellNames(CCLEsmall)




cleanEx()
nameEx("checkPSetStructure")
### * checkPSetStructure

flush(stderr()); flush(stdout())

### Name: checkPSetStructure
### Title: A function to verify the structure of a PharmacoSet
### Aliases: checkPSetStructure

### ** Examples

data(CCLEsmall)

checkPSetStructure(CCLEsmall)




cleanEx()
nameEx("computeABC")
### * computeABC

flush(stderr()); flush(stdout())

### Name: computeABC
### Title: Fits dose-response curves to data given by the user and returns
###   the ABC of the fitted curves.
### Aliases: computeABC

### ** Examples

dose1 <- c("0.0025","0.008","0.025","0.08","0.25","0.8","2.53","8") 
viability1 <- c("108.67","111","102.16","100.27","90","87","74","57")
dose2 <- c("0.0025","0.008","0.025","0.08","0.25","0.8","2.53","8") 
viability2 <- c("100.94","112.5","86","104.16","75","68","48","29")
computeABC(dose1, dose2, viability1, viability2)




cleanEx()
nameEx("computeAUC")
### * computeAUC

flush(stderr()); flush(stdout())

### Name: computeAUC
### Title: Computes the AUC for a Drug Dose Viability Curve
### Aliases: computeAUC

### ** Examples

dose <- c("0.0025","0.008","0.025","0.08","0.25","0.8","2.53","8") 
viability <- c("108.67","111","102.16","100.27","90","87","74","57")
computeAUC(dose, viability)





cleanEx()
nameEx("computeAmax")
### * computeAmax

flush(stderr()); flush(stdout())

### Name: computeAmax
### Title: Fits dose-response curves to data given by the user and returns
###   the Amax of the fitted curve. Amax: 100 - viability at maximum
###   concentarion (in fitted curve)
### Aliases: computeAmax

### ** Examples

dose <- c("0.0025","0.008","0.025","0.08","0.25","0.8","2.53","8") 
viability <- c("108.67","111","102.16","100.27","90","87","74","57")
computeAmax(dose, viability)




cleanEx()
nameEx("computeICn")
### * computeICn

flush(stderr()); flush(stdout())

### Name: computeIC50
### Title: Computes the ICn for any n in 0-100 for a Drug Dose Viability
###   Curve
### Aliases: computeIC50 computeICn

### ** Examples

dose <- c("0.0025","0.008","0.025","0.08","0.25","0.8","2.53","8") 
viability <- c("108.67","111","102.16","100.27","90","87","74","57")
computeIC50(dose, viability)
computeICn(dose, viability, n=10)




cleanEx()
nameEx("computeSlope")
### * computeSlope

flush(stderr()); flush(stdout())

### Name: computeSlope
### Title: Return Slope (normalized slope of the drug response curve) for
###   an experiment of a pSet by taking its concentration and viability as
###   input.
### Aliases: computeSlope

### ** Examples

dose <- c("0.0025","0.008","0.025","0.08","0.25","0.8","2.53","8") 
viability <- c("108.67","111","102.16","100.27","90","87","74","57")
computeSlope(dose, viability)




cleanEx()
nameEx("connectivityScore")
### * connectivityScore

flush(stderr()); flush(stdout())

### Name: connectivityScore
### Title: Function computing connectivity scores between two signatures
### Aliases: connectivityScore

### ** Examples

xValue <- c(1,5,23,4,8,9,2,19,11,12,13)
xSig <- c(0.01, 0.001, .97, 0.01,0.01,0.28,0.7,0.01,0.01,0.01,0.01)
yValue <- c(1,5,10,4,8,19,22,19,11,12,13)
ySig <- c(0.01, 0.001, .97,0.01, 0.01,0.78,0.9,0.01,0.01,0.01,0.01)
xx <- cbind(xValue, xSig)
yy <- cbind(yValue, ySig)
rownames(xx) <- rownames(yy) <- c('1','2','3','4','5','6','7','8','9','10','11')
data.cor <- connectivityScore(xx,yy,method="gwc", gwc.method="spearman", nperm=300)




cleanEx()
nameEx("cosinePerm")
### * cosinePerm

flush(stderr()); flush(stdout())

### Name: cosinePerm
### Title: Computes the cosine similarity and significance using
###   permutation test
### Aliases: cosinePerm

### ** Examples

x <- factor(c(1,2,1,2,1))
y <- factor(c(2,2,1,1,1))
cosinePerm(x, y)




cleanEx()
nameEx("dateCreated")
### * dateCreated

flush(stderr()); flush(stdout())

### Name: dateCreated
### Title: dateCreated Generic
### Aliases: dateCreated

### ** Examples

data(CCLEsmall)
dateCreated(CCLEsmall)




cleanEx()
nameEx("downloadPSet")
### * downloadPSet

flush(stderr()); flush(stdout())

### Name: downloadPSet
### Title: Download a PharmacoSet object
### Aliases: downloadPSet

### ** Examples

if (interactive()){
downloadPSet("CMAP")
}




cleanEx()
nameEx("downloadPertSig")
### * downloadPertSig

flush(stderr()); flush(stdout())

### Name: downloadPertSig
### Title: Download Drug Perturbation Signatures
### Aliases: downloadPertSig

### ** Examples

if (interactive()){
downloadPertSig("CMAP")
}
 



cleanEx()
nameEx("drugDoseResponseCurve")
### * drugDoseResponseCurve

flush(stderr()); flush(stdout())

### Name: drugDoseResponseCurve
### Title: Plot drug response curve of a given drug and a given cell for a
###   list of pSets (objects of the PharmacoSet class).
### Aliases: drugDoseResponseCurve

### ** Examples

if (interactive()) {
drugDoseResponseCurve(concentrations=list("Experiment 1"=c(.008, .04, .2, 1)),
 viabilities=list(c(100,50,30,1)), plot.type="Both")
}




cleanEx()
nameEx("drugInfo-set")
### * drugInfo-set

flush(stderr()); flush(stdout())

### Name: drugInfo<-
### Title: drugInfo<- Generic
### Aliases: drugInfo<-

### ** Examples

data(CCLEsmall)
drugInfo(CCLEsmall) <- drugInfo(CCLEsmall)




cleanEx()
nameEx("drugInfo")
### * drugInfo

flush(stderr()); flush(stdout())

### Name: drugInfo
### Title: drugInfo Generic
### Aliases: drugInfo

### ** Examples

data(CCLEsmall)
drugInfo(CCLEsmall)




cleanEx()
nameEx("drugNames-set")
### * drugNames-set

flush(stderr()); flush(stdout())

### Name: drugNames<-
### Title: drugNames<- Generic
### Aliases: drugNames<-

### ** Examples

data(CCLEsmall)
drugNames(CCLEsmall) <- drugNames(CCLEsmall)




cleanEx()
nameEx("drugNames")
### * drugNames

flush(stderr()); flush(stdout())

### Name: drugNames
### Title: drugNames Generic
### Aliases: drugNames

### ** Examples

data(CCLEsmall)
drugNames(CCLEsmall)




cleanEx()
nameEx("drugPerturbationSig")
### * drugPerturbationSig

flush(stderr()); flush(stdout())

### Name: drugPerturbationSig
### Title: Creates a signature representing gene expression (or other
###   molecular profile) change induced by administrating a drug, for use
###   in drug effect analysis.
### Aliases: drugPerturbationSig

### ** Examples

data(CMAPsmall)
drug.perturbation <- drugPerturbationSig(CMAPsmall, mDataType="rna", nthread=1)
print(drug.perturbation)




cleanEx()
nameEx("drugSensitivitySig")
### * drugSensitivitySig

flush(stderr()); flush(stdout())

### Name: drugSensitivitySig
### Title: Creates a signature representing the association between gene
###   expression (or other molecular profile) and drug dose response, for
###   use in drug sensitivity analysis.
### Aliases: drugSensitivitySig

### ** Examples

data(GDSCsmall)
drug.sensitivity <- drugSensitivitySig(GDSCsmall, mDataType="rna", 
             nthread=1, features = fNames(GDSCsmall, "rna")[1])
print(drug.sensitivity)




cleanEx()
nameEx("fNames")
### * fNames

flush(stderr()); flush(stdout())

### Name: fNames
### Title: fNames Generic
### Aliases: fNames

### ** Examples

data(CCLEsmall)
fNames(CCLEsmall, "rna")




cleanEx()
nameEx("featureInfo-set")
### * featureInfo-set

flush(stderr()); flush(stdout())

### Name: featureInfo<-
### Title: featureInfo<- Generic
### Aliases: featureInfo<-

### ** Examples

data(CCLEsmall)
featureInfo(CCLEsmall, "rna") <- featureInfo(CCLEsmall, "rna")




cleanEx()
nameEx("featureInfo")
### * featureInfo

flush(stderr()); flush(stdout())

### Name: featureInfo
### Title: featureInfo Generic
### Aliases: featureInfo

### ** Examples

data(CCLEsmall)
featureInfo(CCLEsmall, "rna")




cleanEx()
nameEx("filterNoisyCurves")
### * filterNoisyCurves

flush(stderr()); flush(stdout())

### Name: filterNoisyCurves
### Title: Viability measurements in dose-reponse curves must remain stable
###   or decrease monotonically reflecting response to the drug being
###   tested. filterNoisyCurves flags dose-response curves that strongly
###   violate these assumptions.
### Aliases: filterNoisyCurves

### ** Examples

data(GDSCsmall)
filterNoisyCurves(GDSCsmall)




cleanEx()
nameEx("gwc")
### * gwc

flush(stderr()); flush(stdout())

### Name: gwc
### Title: Calculate the gwc score between two vectors, using either a
###   weighted spearman or pearson correlation
### Aliases: gwc

### ** Examples

data(CCLEsmall)
x <- molecularProfiles(CCLEsmall,"rna")[,1]
y <- molecularProfiles(CCLEsmall,"rna")[,2]
x_p <- rep(0.05, times=length(x))
y_p <- rep(0.05, times=length(y))
names(x_p) <- names(x)
names(y_p) <- names(y)
gwc(x,x_p,y,y_p, nperm=100)




cleanEx()
nameEx("intersectList")
### * intersectList

flush(stderr()); flush(stdout())

### Name: intersectList
### Title: Utility to find the intersection between a list of more than two
###   vectors or lists
### Aliases: intersectList

### ** Examples

list1 <- list('a', 'b', 'c')
list2 <- list('a', 'c')
list3 <- list('a', 'c', 'd')
listAll <- intersectList(list1, list2, list3)
listAll




cleanEx()
nameEx("intersectPSet")
### * intersectPSet

flush(stderr()); flush(stdout())

### Name: intersectPSet
### Title: Intersects objects of the PharmacoSet class, subsetting them to
###   the common drugs and/or cell lines as selected by the user.
### Aliases: intersectPSet

### ** Examples

data(GDSCsmall)
data(CCLEsmall)
common <- intersectPSet(list('GDSC'=GDSCsmall,
  'CCLE'=CCLEsmall), intersectOn = c("drugs", "cell.lines"))
common$CGP
common$CCLE





cleanEx()
nameEx("logLogisticRegression")
### * logLogisticRegression

flush(stderr()); flush(stdout())

### Name: logLogisticRegression
### Title: Fits curves of the form E = E_inf + (1 - E_inf)/(1 +
###   (c/EC50)^HS) to dose-response data points (c, E) given by the user
###   and returns a vector containing estimates for HS, E_inf, and EC50.
### Aliases: logLogisticRegression

### ** Examples

dose <- c("0.0025","0.008","0.025","0.08","0.25","0.8","2.53","8") 
viability <- c("108.67","111","102.16","100.27","90","87","74","57")
computeAUC(dose, viability)





cleanEx()
nameEx("mDataNames")
### * mDataNames

flush(stderr()); flush(stdout())

### Name: mDataNames
### Title: mDataNames
### Aliases: mDataNames

### ** Examples

data(CCLEsmall)
mDataNames(CCLEsmall)




cleanEx()
nameEx("mcc")
### * mcc

flush(stderr()); flush(stdout())

### Name: mcc
### Title: Compute a Mathews Correlation Coefficient
### Aliases: mcc

### ** Examples

x <- factor(c(1,2,1,2,3,1))
y <- factor(c(2,1,1,1,2,2))
mcc(x,y)




cleanEx()
nameEx("molecularProfiles-set")
### * molecularProfiles-set

flush(stderr()); flush(stdout())

### Name: molecularProfiles<-
### Title: molecularProfiles<- Generic
### Aliases: molecularProfiles<-

### ** Examples

data(CCLEsmall)
molecularProfiles(CCLEsmall, "rna") <- molecularProfiles(CCLEsmall, "rna")




cleanEx()
nameEx("molecularProfiles")
### * molecularProfiles

flush(stderr()); flush(stdout())

### Name: molecularProfiles
### Title: molecularProfiles Generic
### Aliases: molecularProfiles

### ** Examples

data(CCLEsmall)
molecularProfiles(CCLEsmall, "rna")




cleanEx()
nameEx("pSetName")
### * pSetName

flush(stderr()); flush(stdout())

### Name: pSetName
### Title: pSetName Generic
### Aliases: pSetName

### ** Examples

data(CCLEsmall)
pSetName(CCLEsmall)




cleanEx()
nameEx("pertNumber-set")
### * pertNumber-set

flush(stderr()); flush(stdout())

### Name: pertNumber<-
### Title: pertNumber<- Generic
### Aliases: pertNumber<-

### ** Examples

data(CCLEsmall)
pertNumber(CCLEsmall) <- pertNumber(CCLEsmall)




cleanEx()
nameEx("pertNumber")
### * pertNumber

flush(stderr()); flush(stdout())

### Name: pertNumber
### Title: pertNumber Generic
### Aliases: pertNumber

### ** Examples

data(CCLEsmall)
pertNumber(CCLEsmall)




cleanEx()
nameEx("phenoInfo-set")
### * phenoInfo-set

flush(stderr()); flush(stdout())

### Name: phenoInfo<-
### Title: phenoInfo<- Generic
### Aliases: phenoInfo<-

### ** Examples


data(CCLEsmall)
phenoInfo(CCLEsmall, mDataType="rna") <- phenoInfo(CCLEsmall, mDataType="rna")




cleanEx()
nameEx("phenoInfo")
### * phenoInfo

flush(stderr()); flush(stdout())

### Name: phenoInfo
### Title: phenoInfo Generic
### Aliases: phenoInfo

### ** Examples

data(CCLEsmall)
phenoInfo(CCLEsmall, mDataType="rna")




cleanEx()
nameEx("sensNumber-set")
### * sensNumber-set

flush(stderr()); flush(stdout())

### Name: sensNumber<-
### Title: sensNumber<- Generic
### Aliases: sensNumber<-

### ** Examples

data(CCLEsmall)
sensNumber(CCLEsmall) <- sensNumber(CCLEsmall)




cleanEx()
nameEx("sensNumber")
### * sensNumber

flush(stderr()); flush(stdout())

### Name: sensNumber
### Title: sensNumber Generic
### Aliases: sensNumber

### ** Examples

data(CCLEsmall)
sensNumber(CCLEsmall)




cleanEx()
nameEx("sensitivityInfo-set")
### * sensitivityInfo-set

flush(stderr()); flush(stdout())

### Name: sensitivityInfo<-
### Title: sensitivityInfo<- Generic
### Aliases: sensitivityInfo<-

### ** Examples

data(CCLEsmall)
sensitivityInfo(CCLEsmall) <- sensitivityInfo(CCLEsmall)




cleanEx()
nameEx("sensitivityInfo")
### * sensitivityInfo

flush(stderr()); flush(stdout())

### Name: sensitivityInfo
### Title: sensitivityInfo Generic
### Aliases: sensitivityInfo

### ** Examples

data(CCLEsmall)
sensitivityInfo(CCLEsmall)




cleanEx()
nameEx("sensitivityMeasures")
### * sensitivityMeasures

flush(stderr()); flush(stdout())

### Name: sensitivityMeasures
### Title: sensitivityMeasures Generic
### Aliases: sensitivityMeasures

### ** Examples

data(CCLEsmall)
sensitivityMeasures(CCLEsmall)




cleanEx()
nameEx("sensitivityProfiles-set")
### * sensitivityProfiles-set

flush(stderr()); flush(stdout())

### Name: sensitivityProfiles<-
### Title: sensitivityProfiles<- Generic
### Aliases: sensitivityProfiles<-

### ** Examples

data(CCLEsmall)
sensitivityProfiles(CCLEsmall) <- sensitivityProfiles(CCLEsmall)




cleanEx()
nameEx("sensitivityProfiles")
### * sensitivityProfiles

flush(stderr()); flush(stdout())

### Name: sensitivityProfiles
### Title: sensitivityProfiles Generic
### Aliases: sensitivityProfiles

### ** Examples

data(CCLEsmall)
sensitivityProfiles(CCLEsmall)




cleanEx()
nameEx("show-PharmacoSet-method")
### * show-PharmacoSet-method

flush(stderr()); flush(stdout())

### Name: show,PharmacoSet-method
### Title: Show a PharamcoSet
### Aliases: show,PharmacoSet-method

### ** Examples

data(CCLEsmall)
CCLEsmall




cleanEx()
nameEx("show-PharmacoSig-method")
### * show-PharmacoSig-method

flush(stderr()); flush(stdout())

### Name: show,PharmacoSig-method
### Title: Show PharmacoGx Signatures
### Aliases: show,PharmacoSig-method

### ** Examples

data(GDSCsmall)
drug.sensitivity <- drugSensitivitySig(GDSCsmall, mDataType="rna", 
             nthread=1, features = fNames(GDSCsmall, "rna")[1])
drug.sensitivity




cleanEx()
nameEx("showSigAnnot")
### * showSigAnnot

flush(stderr()); flush(stdout())

### Name: showSigAnnot
### Title: Show the Annotations of a signature object
### Aliases: showSigAnnot

### ** Examples

data(GDSCsmall)
drug.sensitivity <- drugSensitivitySig(GDSCsmall, mDataType="rna", 
             nthread=1, features = fNames(GDSCsmall, "rna")[1])
showSigAnnot(drug.sensitivity)




cleanEx()
nameEx("subsetTo")
### * subsetTo

flush(stderr()); flush(stdout())

### Name: subsetTo
### Title: A function to subset a PharmacoSet to data containing only
###   specified drugs, cells and genes
### Aliases: subsetTo

### ** Examples

data(CCLEsmall)
CCLEdrugs  <- drugNames(CCLEsmall)
CCLEcells <- cellNames(CCLEsmall)
PSet <- subsetTo(CCLEsmall,drugs = CCLEdrugs[1], cells = CCLEcells[1])
PSet




cleanEx()
nameEx("summarizeMolecularProfiles")
### * summarizeMolecularProfiles

flush(stderr()); flush(stdout())

### Name: summarizeMolecularProfiles
### Title: Takes molecular data from a PharmacoSet, and summarises them
###   into one entry per drug
### Aliases: summarizeMolecularProfiles

### ** Examples

data(GDSCsmall)
GDSCsmall <- summarizeMolecularProfiles(GDSCsmall,
                    mDataType = "rna", cell.lines=cellNames(GDSCsmall),
                    summary.stat = 'median', fill.missing = TRUE, verbose=TRUE)
GDSCsmall




cleanEx()
nameEx("summarizeSensitivityProfiles")
### * summarizeSensitivityProfiles

flush(stderr()); flush(stdout())

### Name: summarizeSensitivityProfiles
### Title: Takes the sensitivity data from a PharmacoSet, and summarises
###   them into a drug vs cell line table
### Aliases: summarizeSensitivityProfiles

### ** Examples

data(GDSCsmall)
GDSCauc <- summarizeSensitivityProfiles(GDSCsmall, sensitivity.measure='auc_published')




cleanEx()
nameEx("symSetDiffList")
### * symSetDiffList

flush(stderr()); flush(stdout())

### Name: symSetDiffList
### Title: Utility to find the symmetric set difference of a list of two or
###   more vectors or lists
### Aliases: symSetDiffList

### ** Examples

list1 <- list('a', 'b', 'c')
list2 <- list('a', 'c')
list3 <- list('a', 'c', 'd')
listAll <- symSetDiffList(list1, list2, list3)
listAll




cleanEx()
nameEx("unionList")
### * unionList

flush(stderr()); flush(stdout())

### Name: unionList
### Title: Utility to find the union between a list of more than two
###   vectors or lists
### Aliases: unionList

### ** Examples

list1 <- list('a', 'b')
list2 <- list('a', 'c')
list3 <- list('c', 'd')
listAll <- unionList(list1, list2, list3)
listAll




### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
