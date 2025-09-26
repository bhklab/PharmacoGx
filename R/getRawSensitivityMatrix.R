##TODO:: Add function documentation
getRawSensitivityMatrix <-
  function(pSet, cell.id, drug.id, max.conc, quality) {
    sinfo <- sensitivityInfo(pSet)
    cond_idx <- sinfo$sampleid == cell.id
    if (!missing(quality) && "quality" %in% colnames(sinfo)) {
      cond_idx <- cond_idx & sinfo$quality == quality
    }
    if (!missing(max.conc) && "max.conc" %in% colnames(sinfo)) {
      mc <- if (length(max.conc) > 1) paste(max.conc, collapse = "///") else max.conc
      cond_idx <- cond_idx & sinfo$max.conc == mc
    }
    if (length(drug.id) > 1) {
      drug.id <- paste(drug.id, collapse = "///")
    }
    cond_idx <- cond_idx & sinfo$treatmentid == drug.id

    exp.id <- which(cond_idx)

    sensitivity.raw.matrix <- list()
    if (length(exp.id) > 0) {
      for (i in seq_len(length(exp.id))) {
        if (grepl("///", drug.id, fixed = TRUE)) {
          all.exp.id <- which(
            sensitivityInfo(pSet)$combination.exp.id ==
              sensitivityInfo(pSet)[exp.id[i], "combination.exp.id"]
          )
          drug.1 <- which(
            sensitivityInfo(pSet)[all.exp.id, "treatmentid"] ==
              unlist(strsplit(drug.id, split = "///"))[1]
          )
          drug.2 <- which(
            sensitivityInfo(pSet)[all.exp.id, "treatmentid"] ==
              unlist(strsplit(drug.id, split = "///"))[2]
          )
          drug.1.doses <- length(which(
            !is.na(sensitivityRaw(pSet)[all.exp.id[drug.1], , "Dose"])
          ))
          drug.2.doses <- length(which(
            !is.na(sensitivityRaw(pSet)[all.exp.id[drug.2], , "Dose"])
          ))

          tt <- matrix(NA, ncol = drug.2.doses, nrow = drug.1.doses)
          colnames(tt) <- sensitivityRaw(pSet)[
            all.exp.id[drug.2],
            seq_len(drug.2.doses),
            "Dose"
          ]
          rownames(tt) <- sensitivityRaw(pSet)[
            all.exp.id[drug.1],
            seq_len(drug.1.doses),
            "Dose"
          ]
          tt[, 1] <- sensitivityRaw(pSet)[
            all.exp.id[drug.1],
            seq_len(drug.1.doses),
            "Viability"
          ]
          tt[1, ] <- sensitivityRaw(pSet)[
            all.exp.id[drug.2],
            seq_len(drug.2.doses),
            "Viability"
          ]
          tt[2:nrow(tt), 2:ncol(tt)] <- sensitivityRaw(pSet)[
            exp.id[i],
            ,
            "Viability"
          ]
          sensitivity.raw.matrix[[rownames(sensitivityInfo(pSet))[exp.id[
            i
          ]]]] <- tt
        } else {
          sensitivity.raw.matrix[[rownames(sensitivityInfo(pSet))[exp.id[
            i
          ]]]] <- sensitivityRaw(pSet)[exp.id[i], , ]
        }
      }
    }
    return(sensitivity.raw.matrix)
  }
