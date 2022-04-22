##TODO:: Add function documentation
getRawSensitivityMatrix <-
  function(pSet, cell.id, drug.id, max.conc, quality) {
    cond <- "sensitivityInfo(pSet)$sampleid == cell.id"
    if(!missing(quality)) {
      if(is.element("quality", colnames(sensitivityInfo(pSet)))) {
        cond <- paste(cond, "sensitivityInfo(pSet)$quality == quality", sep=" & ")
      }
    }
    if(!missing(max.conc)) {
        if(is.element("max.conc", colnames(sensitivityInfo(pSet)))) {
          if(length(max.conc) > 1) {
            max.conc <- paste(max.conc, collapse="///")
        }
        cond <- paste(cond, "sensitivityInfo(pSet)$max.conc == max.conc", sep=" & ")
      }
    }
    if(length(drug.id) > 1) {
      drug.id <- paste(drug.id, collapse="///")
    }
    cond <- paste(cond, "sensitivityInfo(pSet)$treatmentid == drug.id", sep=" & ")

    exp.id <- which(eval(parse(text=cond)))

    sensitivity.raw.matrix <- list()
    if(length(exp.id) > 0) {
      for(i in seq_len(length(exp.id))){
        if(length(grep("///", drug.id)) > 0) {
          all.exp.id <- which(sensitivityInfo(pSet)$combination.exp.id == sensitivityInfo(pSet)[exp.id[i], "combination.exp.id"])
          drug.1 <- which(sensitivityInfo(pSet)[all.exp.id, "treatmentid"] == unlist(strsplit(drug.id, split="///"))[1])
          drug.2 <- which(sensitivityInfo(pSet)[all.exp.id, "treatmentid"] == unlist(strsplit(drug.id, split="///"))[2])
          drug.1.doses <- length(which(!is.na(sensitivityRaw(pSet)[all.exp.id[drug.1], , "Dose"])))
          drug.2.doses <- length(which(!is.na(sensitivityRaw(pSet)[all.exp.id[drug.2], , "Dose"])))

          tt <- matrix(NA, ncol=drug.2.doses, nrow=drug.1.doses)
          colnames(tt) <- sensitivityRaw(pSet)[all.exp.id[drug.2], seq_len(drug.2.doses), "Dose"]
          rownames(tt) <- sensitivityRaw(pSet)[all.exp.id[drug.1], seq_len(drug.1.doses), "Dose"]
          tt[ ,1] <- sensitivityRaw(pSet)[all.exp.id[drug.1], seq_len(drug.1.doses), "Viability"]
          tt[1, ] <- sensitivityRaw(pSet)[all.exp.id[drug.2], seq_len(drug.2.doses), "Viability"]
          tt[2:nrow(tt), 2:ncol(tt)] <- sensitivityRaw(pSet)[exp.id[i], , "Viability"]
          sensitivity.raw.matrix[[rownames(sensitivityInfo(pSet))[exp.id[i]]]] <- tt
        }else{
          sensitivity.raw.matrix[[rownames(sensitivityInfo(pSet))[exp.id[i]]]] <- sensitivityRaw(pSet)[exp.id[i], , ]
        }
      }
    }
    return(sensitivity.raw.matrix)
  }
