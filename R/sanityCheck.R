.normalizeHillPars <- function(
  hill_fit,
  conc_as_log = FALSE,
  viability_as_pct = TRUE
) {
  if (is.null(hill_fit)) {
    return(NULL)
  }

  rsq <- attr(hill_fit, "Rsquare", exact = TRUE)
  hill <- unlist(hill_fit, use.names = TRUE)

  if (!length(hill)) {
    stop("`Hill_fit` must contain at least one parameter.")
  }

  lower_names <- tolower(names(hill))
  if (all(is.na(lower_names))) {
    lower_names <- rep("", length(hill))
  }

  get_param <- function(candidates, fallback_index = NA_integer_, default = NA_real_) {
    idx <- match(candidates, lower_names)
    idx <- idx[!is.na(idx)][1]
    if (!is.na(idx)) {
      return(as.numeric(hill[idx]))
    }
    if (!is.na(fallback_index) && fallback_index <= length(hill)) {
      return(as.numeric(hill[fallback_index]))
    }
    default
  }

  hs <- get_param(c("hs", "hill", "hill_slope"), fallback_index = 1)
  e0 <- get_param(c("e0", "etop"), default = ifelse(viability_as_pct, 100, 1))
  einf <- get_param(c("e_inf", "einf", "einfty"), fallback_index = if (length(hill) >= 3) 2 else NA)
  ec50 <- get_param(c("ec50", "logec50", "log10ec50"), fallback_index = if (length(hill) >= 3) 3 else NA)

  if (any(!is.finite(c(hs, e0, einf, ec50)))) {
    stop("Unable to parse `Hill_fit` parameters into a numeric vector.")
  }

  if (!conc_as_log) {
    if (ec50 <= 0) {
      stop("EC50 must be positive when `conc_as_log = FALSE`.")
    }
    log_ec50 <- log10(ec50)
  } else {
    log_ec50 <- ec50
  }

  if (viability_as_pct) {
    e0 <- e0 / 100
    einf <- einf / 100
  }

  pars <- c(
    HS = hs,
    E0 = e0,
    E_inf = einf,
    log10EC50 = log_ec50
  )
  attr(pars, "Rsquare") <- rsq
  pars
}

.normalizeBiphasicPars <- function(
  hill_fit,
  conc_as_log = FALSE,
  viability_as_pct = TRUE
) {
  if (is.null(hill_fit)) {
    return(NULL)
  }

  rsq <- attr(hill_fit, "Rsquare", exact = TRUE)
  hill <- unlist(hill_fit, use.names = TRUE)
  lower_names <- tolower(names(hill))

  get_param <- function(candidates, fallback_index = NA_integer_, default = NA_real_) {
    idx <- match(candidates, lower_names)
    idx <- idx[!is.na(idx)][1]
    if (!is.na(idx)) {
      return(as.numeric(hill[idx]))
    }
    if (!is.na(fallback_index) && fallback_index <= length(hill)) {
      return(as.numeric(hill[fallback_index]))
    }
    default
  }

  hs1 <- get_param(c("hs1"), fallback_index = 1, default = 1)
  e0 <- get_param(c("e0", "etop"), fallback_index = 2, default = ifelse(viability_as_pct, 100, 1))
  einf1 <- get_param(c("e_inf1", "einf1"), fallback_index = 3, default = 0)
  hs2 <- get_param(c("hs2"), fallback_index = 4, default = 1)
  einf2 <- get_param(c("e_inf2", "einf2"), fallback_index = 5, default = 0)
  ec50_1 <- get_param(c("ec50_1", "logec50_1", "log10ec50_1"), fallback_index = 6, default = NA_real_)
  ec50_2 <- get_param(c("ec50_2", "logec50_2", "log10ec50_2"), fallback_index = 7, default = NA_real_)
  frac <- get_param(c("frac"), fallback_index = 8, default = 0.5)

  params <- c(hs1, e0, einf1, hs2, einf2, ec50_1, ec50_2, frac)
  if (any(!is.finite(params))) {
    stop("Unable to parse biphasic parameters from `Hill_fit`.")
  }

  if (!conc_as_log) {
    if (ec50_1 <= 0 || ec50_2 <= 0) {
      stop("EC50 parameters must be positive when `conc_as_log = FALSE`.")
    }
    log_ec50_1 <- log10(ec50_1)
    log_ec50_2 <- log10(ec50_2)
  } else {
    log_ec50_1 <- ec50_1
    log_ec50_2 <- ec50_2
  }

  if (viability_as_pct) {
    e0 <- e0 / 100
    einf1 <- einf1 / 100
    einf2 <- einf2 / 100
  }

  pars <- c(
    HS1 = hs1,
    E0 = e0,
    E_inf1 = einf1,
    HS2 = hs2,
    E_inf2 = einf2,
    log10EC50_1 = log_ec50_1,
    log10EC50_2 = log_ec50_2,
    Frac = frac
  )
  attr(pars, "Rsquare") <- rsq
  pars
}

sanitizeInput <- function(
  conc,
  viability,
  Hill_fit,
  conc_as_log = FALSE,
  viability_as_pct = TRUE,
  trunc = TRUE,
  verbose = TRUE # Set to 2 to see debug printouts
) {
  if (is.logical(conc_as_log) == FALSE) {
    print(conc_as_log)
    stop("'conc_as_log' is not a logical.")
  }

  if (is.logical(viability_as_pct) == FALSE) {
    print(viability_as_pct)
    stop("'viability_as_pct' is not a logical.")
  }

  if (is.logical(trunc) == FALSE) {
    print(trunc)
    stop("'trunc' is not a logical.")
  }
  if (!is.finite(verbose)) {
    stop("'verbose' should be a logical (or numerical) argument.")
  }
  if (!missing(viability) && !missing(conc) && missing(Hill_fit)) {
    if (length(conc) != length(viability)) {
      if (verbose == 2) {
        print(conc)
        print(viability)
      }
      stop(
        "Log concentration vector is not of same length as viability vector."
      )
    }
    if (any(is.na(conc) & (!is.na(viability)))) {
      warning(
        "Missing concentrations with non-missing viability values encountered. Removing viability values corresponding to those concentrations"
      )

      myx <- !is.na(conc)
      conc <- as.numeric(conc[myx])
      viability <- as.numeric(viability[myx])
    }
    if (any((!is.na(conc)) & is.na(viability))) {
      warning(
        "Missing viability with non-missing concentrations values encountered. Removing concentrations values corresponding to those viabilities"
      )
      myx <- !is.na(viability)
      conc <- as.numeric(conc[myx])
      viability <- as.numeric(viability[myx])
    }
    conc <- as.numeric(conc[!is.na(conc)])
    viability <- as.numeric(viability[!is.na(viability)])

    #CHECK THAT FUNCTION INPUTS ARE APPROPRIATE
    if (prod(is.finite(conc)) != 1) {
      print(conc)
      stop("Concentration vector contains elements which are not real numbers.")
    }

    if (prod(is.finite(viability)) != 1) {
      print(viability)
      stop("Viability vector contains elements which are not real numbers.")
    }

    if (min(viability) < 0) {
      if (verbose) {
        warning("Warning: Negative viability data.")
      }
    }

    if (max(viability) > (1 + 99 * viability_as_pct)) {
      if (verbose) {
        warning("Warning: Viability data exceeds negative control.")
      }
    }

    if (conc_as_log == FALSE && min(conc) < 0) {
      if (verbose == 2) {
        print(conc)
        print(conc_as_log)
      }
      stop(
        "Negative concentrations encountered. Concentration data may be inappropriate, or 'conc_as_log' flag may be set incorrectly."
      )
    }

    if (viability_as_pct == TRUE && max(viability) < 5) {
      warning("Warning: 'viability_as_pct' flag may be set incorrectly.")
      if (verbose == 2) {
        print(viability)
        print(viability_as_pct)
      }
    }

    if (viability_as_pct == FALSE && max(viability) > 5) {
      warning("Warning: 'viability_as_pct' flag may be set incorrectly.")
      if (verbose == 2) {
        print(viability)
        print(viability_as_pct)
      }
    }

    if (is.unsorted(conc)) {
      warning(
        "Concentration Values were unsorted. Sorting concentration and ordering viability in same order"
      )
      myx <- order(conc)
      conc <- conc[myx]
      viability <- viability[myx]
    }

    #CONVERT DOSE-RESPONSE DATA TO APPROPRIATE INTERNAL REPRESENTATION
    if (conc_as_log == FALSE) {
      ii <- which(conc == 0)
      if (length(ii) > 0) {
        conc <- conc[-ii]
        viability <- viability[-ii]
      }

      log_conc <- log10(conc)
    } else {
      log_conc <- conc
    }

    if (viability_as_pct == TRUE) {
      viability <- viability / 100
    }
    if (trunc) {
      viability = pmin(as.numeric(viability), 1)
      viability = pmax(as.numeric(viability), 0)
    }

    return(list("log_conc" = log_conc, "viability" = viability))
  }
  if (!missing(Hill_fit) && missing(viability)) {
    raw_fit <- Hill_fit
    lower_names <- tolower(names(unlist(raw_fit, use.names = TRUE)))
    if (any(lower_names %in% c("hs1", "e_inf1", "einf1")) || length(unlist(raw_fit)) >= 8) {
      Hill_fit <- .normalizeBiphasicPars(
        hill_fit = raw_fit,
        conc_as_log = conc_as_log,
        viability_as_pct = viability_as_pct
      )
    } else {
      Hill_fit <- .normalizeHillPars(
        hill_fit = raw_fit,
        conc_as_log = conc_as_log,
        viability_as_pct = viability_as_pct
      )
    }

    if (!viability_as_pct && Hill_fit[["E_inf"]] > 1) {
      warning("'viability_as_pct' flag may be set incorrectly.")
    }

    if (missing(conc)) {
      return(list("Hill_fit" = Hill_fit))
    }

    conc <- as.numeric(conc[!is.na(conc)])

    if (prod(is.finite(conc)) != 1) {
      print(conc)
      stop("Concentration vector contains elements which are not real numbers.")
    }
    if (conc_as_log == FALSE && min(conc) < 0) {
      print(conc)
      print(conc_as_log)
      stop(
        "Negative concentrations encountered. Concentration data may be inappropriate, or 'conc_as_log' flag may be set incorrectly."
      )
    }

    if (conc_as_log == FALSE) {
      ii <- which(conc == 0)
      if (length(ii) > 0) {
        conc <- conc[-ii]
      }
      log_conc <- log10(conc)
    } else {
      log_conc <- conc
    }
    if (is.unsorted(conc)) {
      myx <- order(conc)
      conc <- conc[myx]
    }
    return(list("Hill_fit" = Hill_fit, "log_conc" = log_conc))
  }
  if (!missing(Hill_fit) && !missing(viability)) {
    stop(
      "Please pass in only one of 'Hill_fit' and 'viability', it is unclear which to use in the computation."
    )
  }
  if (missing(Hill_fit) && missing(viability)) {
    stop("Both 'Hill_fit' and 'viability' missing, please pass in some data!")
  }
}
