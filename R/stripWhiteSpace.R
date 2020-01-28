stripWhiteSpace <- function (str, method=c("both", "head", "tail")) {
  method <- match.arg(method)
  str2 <- NULL
  if (length(str) == 1) {
    switch (method,
      "both" = {
       str2 <- gsub("^[ \t]+", "", str)
       str2 <- gsub("[ \t]+$", "", str2)
      },
      "head" = {
        str2 <- gsub("^[ \t]+", "", str)
      },
      "tail" = {
        str2 <- gsub("[ \t]+$", "", str)
      }
    )
    return (str2)
  } else {
    str2 <- sapply(str, stripWhiteSpace, method=method)
    return (str2)
  }
}

## End

