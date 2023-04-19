## ----setup, include=FALSE, cache=FALSE, message = FALSE-----------------------

library("knitr")

### Chunk options: see http://yihui.name/knitr/options/ ###

## Text results
opts_chunk$set(echo = TRUE, warning = TRUE, message = TRUE, include = TRUE)

## Code decoration
opts_chunk$set(tidy = FALSE, comment = NA, highlight = TRUE)


## ----knitcitation, include=FALSE----------------------------------------------
library(knitcitations)
cleanbib()
cite_options(citation_format = "pandoc")

## ----eval=FALSE---------------------------------------------------------------
#  SE@annotation <- "rna"

## ----eval=FALSE---------------------------------------------------------------
#  PharmacoSet(name,
#              molecularProfiles=list(),
#              sample=data.frame(),
#              treatment=data.frame(),
#              sensitivityInfo=data.frame(),
#              sensitivityRaw=array(dim=c(0,0,0)),
#              sensitivityProfiles=matrix(),
#              curationTreatment=data.frame(),
#              curationSample=data.frame(),
#              curationTissue=data.frame(),
#              datasetType=c("sensitivity", "perturbation", "both"),
#              verify = TRUE)

