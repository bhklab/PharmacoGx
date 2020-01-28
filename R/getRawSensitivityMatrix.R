getRawSensitivityMatrix <-
  function(pSet, cell.id, drug.id, max.conc, quality) {
    cond <- "pSet@sensitivity$info$cellid == cell.id"
    if(!missing(quality)) {
      if(is.element("quality", colnames(pSet@sensitivity$info))) {
        cond <- paste(cond, "pSet@sensitivity$info$quality == quality", sep=" & ")
      }
    }
    if(!missing(max.conc)) {
        if(is.element("max.conc", colnames(pSet@sensitivity$info))) {
          if(length(max.conc) > 1) {
            max.conc <- paste(max.conc, collapse="///")
        }
        cond <- paste(cond, "pSet@sensitivity$info$max.conc == max.conc", sep=" & ")
      }
    }
    if(length(drug.id) > 1) {
      drug.id <- paste(drug.id, collapse="///")
    }
    cond <- paste(cond, "pSet@sensitivity$info$drugid == drug.id", sep=" & ")
    
    exp.id <- which(eval(parse(text=cond)))
    
    sensitivity.raw.matrix <- list()
    if(length(exp.id) > 0) {
      for(i in 1:length(exp.id)){
        if(length(grep("///", drug.id)) > 0) {
          all.exp.id <- which(pSet@sensitivity$info$combination.exp.id == pSet@sensitivity$info[exp.id[i], "combination.exp.id"])
          drug.1 <- which(pSet@sensitivity$info[all.exp.id, "drugid"] == unlist(strsplit(drug.id, split="///"))[1])
          drug.2 <- which(pSet@sensitivity$info[all.exp.id, "drugid"] == unlist(strsplit(drug.id, split="///"))[2])
          drug.1.doses <- length(which(!is.na(pSet@sensitivity$raw[all.exp.id[drug.1], , "Dose"])))
          drug.2.doses <- length(which(!is.na(pSet@sensitivity$raw[all.exp.id[drug.2], , "Dose"])))
          
          tt <- matrix(NA, ncol=drug.2.doses, nrow=drug.1.doses)
          colnames(tt) <- pSet@sensitivity$raw[all.exp.id[drug.2], 1:drug.2.doses, "Dose"]
          rownames(tt) <- pSet@sensitivity$raw[all.exp.id[drug.1], 1:drug.1.doses, "Dose"]
          tt[ ,1] <- pSet@sensitivity$raw[all.exp.id[drug.1], 1:drug.1.doses, "Viability"]
          tt[1, ] <- pSet@sensitivity$raw[all.exp.id[drug.2], 1:drug.2.doses, "Viability"]
          tt[2:nrow(tt), 2:ncol(tt)] <- pSet@sensitivity$raw[exp.id[i], , "Viability"]
          sensitivity.raw.matrix[[rownames(pSet@sensitivity$info)[exp.id[i]]]] <- tt
        }else{
          sensitivity.raw.matrix[[rownames(pSet@sensitivity$info)[exp.id[i]]]] <- pSet@sensitivity$raw[exp.id[i], , ]
        }
      }
    }
    return(sensitivity.raw.matrix)
  }

