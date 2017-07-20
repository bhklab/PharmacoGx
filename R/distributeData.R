.distributeData <- function(data, nslaves, split, byRows = TRUE) {
  xx <- list()
  xx[[1]] <- NA #to ensure that the master node does not get any data
  
  if (missing(split)) {
    if (is.matrix(data) || is.data.frame(data) || is.vector(data) || is.list(data)) {
      split <- TRUE
    } else {
      split <- FALSE
    }
  }
  
  if (split) {
    if (is.matrix(data) || is.data.frame(data)) {
      if (byRows) {
        if(nrow(data) < nslaves){
          stop(paste0("Not enough data for all slaves to be kept busy."))
        }
        rowsPerSlave <- floor(nrow(data) / nslaves)
        extraRows <- nrow(data) %% nslaves
        rowIndex <- 1
        for (i in seq_len(extraRows)) {
          xx[[i + 1]] <- data[rowIndex:(rowIndex + rowsPerSlave), , drop=FALSE]
          rowIndex <- rowIndex + rowsPerSlave + 1
        }
        for (i in (extraRows + 1):nslaves) {
          xx[[i + 1]] <- data[rowIndex:(rowIndex + rowsPerSlave - 1), , drop=FALSE]
          rowIndex <- rowIndex + rowsPerSlave
        }
      } else {
        if(ncol(data) < nslaves){
          stop(paste0("Not enough data for all slaves to be kept busy."))
        }
        colsPerSlave <- floor(ncol(data) / nslaves)
        extraCols <- ncol(data) %% nslaves
        colIndex <- 1
        for (i in seq_len(extraCols)) {
          xx[[i + 1]] <- data[, colIndex:(colIndex + colsPerSlave), drop=FALSE]
          colIndex <- colIndex + colsPerSlave + 1
        }
        for (i in (extraCols + 1):nslaves) {
          xx[[i + 1]] <- data[, colIndex:(colIndex + colsPerSlave - 1), drop=FALSE]
          colIndex <- colIndex + colsPerSlave
        }
      }
    } else if (is.vector(data) || is.list(data)) {
      if(length(data) < nslaves){
        stop(paste0("Not enough data for all slaves to be kept busy."))
      }
      elementsPerSlave <- floor(length(data) / nslaves)
      extraElements <- length(data) %% nslaves
      elementIndex <- 1
      for (i in seq_len(extraElements)) {
        xx[[i + 1]] <- data[[elementIndex:(elementIndex + elementsPerSlave)]]
        elementIndex <- elementIndex + elementsPerSlave + 1
      }
      for (i in (extraElements + 1):nslaves) {
        xx[[i + 1]] <- data[[elementIndex:(elementIndex + elementsPerSlave - 1)]]
        elementIndex <- elementIndex + elementsPerSlave
      }
    } else {
      stop("Splitting this type of object is not yet supported by .distributeData.")
    }
  } else {
    for (i in seq_len(nslaves)) {
      xx[[i + 1]] <- data
    }
  }
  
  return(xx)
}