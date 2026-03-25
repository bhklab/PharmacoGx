#' Plot drug response curve of a given drug and a given cell for a list of pSets (objects of the PharmacoSet class).
#'
#' Given a list of PharmacoSets, the function will plot the drug_response curve,
#' for a given drug/cell pair. The y axis of the plot is the viability percentage
#' and x axis is the log transformed concentrations. If more than one pSet is
#' provided, a light gray area would show the common concentration range between pSets.
#' User can ask for type of sensitivity measurment to be shown in the plot legend.
#' The user can also provide a list of their own concentrations and viability values,
#' as in the examples below, and it will be treated as experiments equivalent to values coming
#' from a pset. The names of the concentration list determine the legend labels.
#'
#' @examples
##TODO:: How do you pass PSets to this?
#' if (interactive()) {
#' # Manually enter the plot parameters
#' drugDoseResponseCurve(concentrations=list("Experiment 1"=c(.008, .04, .2, 1)),
#'  viabilities=list(c(100,50,30,1)), plot.type="Both")
#'
#' # Generate a plot from one or more PSets
#' data(GDSCsmall)
#' drugDoseResponseCurve(drug="Doxorubicin", cellline="22RV1", pSets=GDSCsmall)
#' }
#'
#' @param drug `character(1)` A drug name for which the drug response curve should be
#' plotted. If the plot is desirable for more than one pharmaco set, A unique drug id
#' should be provided.
#' @param cellline `character(1)` A cell line name for which the drug response curve should be
#' plotted. If the plot is desirable for more than one pharmaco set, A unique cell id
#' should be provided.
#' @param pSets `list` a list of PharmacoSet objects, for which the function
#' should plot the curves.
#' @param concentrations,viabilities `list` A list of concentrations and viabilities to plot, the function assumes that
#' `concentrations[[i]]` is plotted against `viabilities[[i]]`. The names of the concentration list are used to create the legend labels
#' @param conc_as_log `logical`, if true, assumes that log10-concentration data has been given rather than concentration data,
#' and that log10(ICn) should be returned instead of ICn. Applies only to the concentrations parameter.
#' @param viability_as_pct `logical`, if false, assumes that viability is given as a decimal rather
#' than a percentage, and that E_inf passed in as decimal. Applies only to the viabilities parameter.
#' @param legends.label `numeric` A vector of sensitivity measurment types which could
#' be any combination of  ic50_published, auc_published, auc_recomputed and auc_recomputed_star.
#' A legend will be displayed on the top right of the plot which each line of the legend is
#' the values of requested sensitivity measerments for one of the requested pSets.
#' If this parameter is missed no legend would be provided for the plot.
#' @param ylim `numeric` A vector of two numerical values to be used as ylim of the plot.
#' If this parameter would be missed c(0,100) would be used as the ylim of the plot.
#' @param xlim `numeric` A vector of two numerical values to be used as xlim of the plot.
#' If this parameter would be missed the minimum and maximum comncentrations between all
#' the pSets would be used as plot xlim.
#' @param mycol `numeric` A vector with the same lenght of the pSets parameter which
#' will determine the color of the curve for the pharmaco sets. If this parameter is
#' missed default colors from Rcolorbrewer package will be used as curves color.
#' @param plot.type `character` Plot type which can be the actual one ("Actual") or
#' the one fitted by logl logistic regression ("Fitted") or both of them ("Both").
#' If this parameter is missed by default actual curve is plotted.
#' @param summarize.replicates `logical(1)` If `TRUE`, replicate measurements at
#' identical concentrations are summarized by median viability before plotting
#' and fitting. If replicate experiments have non-identical dose grids, only
#' exact duplicate doses are summarized and a warning is emitted. If `FALSE`,
#' replicate experiments are plotted individually.
#' @param title `character` The title of the graph. If no title is provided, then it defaults to
#' 'Drug':'Cell Line'.
#' @param lwd `numeric` The line width to plot with
#' @param cex `numeric` The cex parameter passed to plot
#' @param cex.main `numeric` The cex.main parameter passed to plot, controls the size of the titles
#' @param legend.loc And argument passable to xy.coords for the position to place the legend.
#' @param trunc `logical(1)` Should the viability values be truncated to lie in \[0-100\] before doing the fitting
#' @param verbose `logical(1)` Should warning messages about the data passed in be printed?
#' @param sample_col `character(1)` The name of the column in the profiles assay that contains the sample IDs.
#' @param treatment_col `character(1)` The name of the column in the profiles assay that contains the treatment IDs.
#'
#' @return Plots to the active graphics device and returns an invisible NULL.
#'
#' @import RColorBrewer
#'
#' @importFrom graphics plot rect points lines legend
#' @importFrom grDevices rgb
# # ' @importFrom magicaxis magaxis
#' @importFrom CoreGx .getSupportVec
#'
#' @export
drugDoseResponseCurve <-
  function(
    drug,
    cellline,
    pSets = list(),
    concentrations = list(),
    viabilities = list(),
    conc_as_log = FALSE,
    viability_as_pct = TRUE,
    trunc = TRUE,
    legends.label = c(
      "ic50_published",
      "gi50_published",
      "auc_published",
      "auc_recomputed",
      "ic50_recomputed"
    ),
    ylim = c(0, 100),
    xlim,
    mycol,
    title,
    plot.type = c("Actual", "Fitted", "Both"),
    summarize.replicates = TRUE,
    lwd = 0.5,
    cex = 0.7,
    cex.main = 0.9,
    legend.loc = "topright",
    verbose = TRUE,
    sample_col = "sampleid",
    treatment_col = "treatmentid"
  ) {
    if (!missing(pSets)) {
      if (!is(pSets, "list")) {
        if (is(pSets, "PharmacoSet")) {
          temp <- name(pSets)
          pSets <- list(pSets)
          names(pSets) <- temp
        } else {
          stop(
            "Type of pSets parameter should be either a pSet or a list of pSets."
          )
        }
      }
    }
    if (!missing(pSets) && (missing(drug) || missing(cellline))) {
      stop("If you pass in a pSet then drug and cellline must be set")
    }
    # } else {
    #   if(missing(drug)){
    #   drug <- "Drug"}
    #   if(missing(cellline))
    #   cellline <- "Cell Line"
    # }
    if (!missing(concentrations)) {
      if (missing(viabilities)) {
        stop("Please pass in the viabilities to Plot with the concentrations.")
      }
      if (!is(concentrations, "list")) {
        if (mode(concentrations) == "numeric") {
          if (mode(viabilities) != "numeric") {
            stop(
              "Passed in 1 vector of concentrations but the viabilities are not numeric!"
            )
          }
          cleanData <- sanitizeInput(
            concentrations,
            viabilities,
            conc_as_log = conc_as_log,
            viability_as_pct = viability_as_pct,
            trunc = trunc,
            verbose = verbose
          )
          concentrations <- 10^cleanData[["log_conc"]]
          concentrations <- list(concentrations)
          viabilities <- 100 * cleanData[["viability"]]
          viabilities <- list(viabilities)
          names(concentrations) <- "Exp1"
          names(viabilities) <- "Exp1"
        } else {
          stop(
            "Mode of concentrations parameter should be either numeric or a list of numeric vectors"
          )
        }
      } else {
        if (length(viabilities) != length(concentrations)) {
          stop(
            "The number of concentration and viability vectors passed in differs"
          )
        }
        if (is.null(names(concentrations))) {
          names(concentrations) <- paste("Exp", seq_len(length(concentrations)))
        }
        for (i in seq_len(length(concentrations))) {
          if (mode(concentrations[[i]]) == "numeric") {
            if (mode(viabilities[[i]]) != "numeric") {
              stop(sprintf(
                "concentrations[[%d]] are numeric but the viabilities[[%d]] are not numeric!",
                i,
                i
              ))
            }
            cleanData <- sanitizeInput(
              concentrations[[i]],
              viabilities[[i]],
              conc_as_log = conc_as_log,
              viability_as_pct = viability_as_pct,
              trunc = trunc,
              verbose = verbose
            )
            concentrations[[i]] <- 10^cleanData[["log_conc"]]
            viabilities[[i]] <- 100 * cleanData[["viability"]]
          } else {
            stop(sprintf(
              "Mode of concentrations[[%d]] parameter should be numeric",
              i
            ))
          }
        }
      }
    }

    plot.type <- match.arg(plot.type)

    .extract_drug_responses <- function(pSet, experiments) {
      response.list <- lapply(experiments, function(exp) {
        drug.responses <- data.frame(
          Dose = as.numeric(as.vector(sensitivityRaw(pSet)[exp, , "Dose"])),
          Viability = as.numeric(as.vector(sensitivityRaw(pSet)[
            exp,
            ,
            "Viability"
          ])),
          stringsAsFactors = FALSE
        )
        valid.idx <- complete.cases(drug.responses) &
          is.finite(drug.responses$Dose) &
          is.finite(drug.responses$Viability) &
          drug.responses$Dose > 0
        drug.responses[valid.idx, , drop = FALSE]
      })
      names(response.list) <- rownames(sensitivityInfo(pSet))[experiments]
      response.list
    }

    .dose_grids_identical <- function(response.list) {
      if (length(response.list) <= 1) {
        return(TRUE)
      }

      reference.dose <- sort(response.list[[1]]$Dose)
      all(vapply(
        response.list[-1],
        function(drug.responses) {
          identical(sort(drug.responses$Dose), reference.dose)
        },
        logical(1)
      ))
    }

    .summarize_responses <- function(response.list) {
      combined.responses <- do.call(rbind, response.list)
      if (is.null(combined.responses) || nrow(combined.responses) == 0) {
        return(data.frame(Dose = numeric(), Viability = numeric()))
      }

      summarized.responses <- stats::aggregate(
        Viability ~ Dose,
        data = combined.responses,
        FUN = function(x) median(as.numeric(x), na.rm = TRUE)
      )
      summarized.responses[order(summarized.responses$Dose), , drop = FALSE]
    }

    for (i in seq_len(length(pSets))) {
      if (is(treatmentResponse(pSets[[i]]), "LongTable")) {
        pSets[[i]] <- subsetByTreatment(pSets[[i]], treatments = drug)
      }
      pSets[[i]] <- subsetBySample(pSets[[i]], samples = cellline)
    }

    doses <- list()
    responses <- list()
    legend.values <- list()
    j <- 0
    pSetNames <- list()
    if (!missing(pSets)) {
      for (i in seq_len(length(pSets))) {
        exp_i <- which(
          sensitivityInfo(pSets[[i]])[, sample_col] == cellline &
            sensitivityInfo(pSets[[i]])[, treatment_col] == drug
        )
        if (length(exp_i) > 0) {
          if (summarize.replicates) {
            pSetNames[[i]] <- name(pSets[[i]])
            response.list <- .extract_drug_responses(pSets[[i]], exp_i)
            identical.dose.grids <- .dose_grids_identical(response.list)
            drug.responses <- .summarize_responses(response.list)

            if (length(response.list) > 1 && !identical.dose.grids) {
              warning(
                sprintf(
                  paste(
                    "Replicate dose grids differ for %s:%s in %s;",
                    "only exact duplicate doses were summarized."
                  ),
                  drug,
                  cellline,
                  name(pSets[[i]])
                )
              )
            }

            doses[[i]] <- drug.responses$Dose
            responses[[i]] <- drug.responses$Viability
            names(doses[[i]]) <- names(responses[[i]]) <- seq_len(length(doses[[
              i
            ]]))
            if (!missing(legends.label)) {
              if (length(legends.label) > 1) {
                legend.values[[i]] <- paste(
                  unlist(lapply(legends.label, function(x) {
                    sprintf(
                      "%s = %s",
                      x,
                      round(
                        as.numeric(sensitivityProfiles(pSets[[i]])[exp_i, x]),
                        digits = 2
                      )
                    )
                  })),
                  collapse = ", "
                )
              } else {
                legend.values[[i]] <- sprintf(
                  "%s = %s",
                  legends.label,
                  round(
                    as.numeric(sensitivityProfiles(pSets[[i]])[
                      exp_i,
                      legends.label
                    ]),
                    digits = 2
                  )
                )
              }
            } else {
              legend.values[i] <- ""
            }
          } else {
            for (exp in exp_i) {
              j <- j + 1
              pSetNames[[j]] <- name(pSets[[i]])

              drug.responses <- as.data.frame(
                cbind(
                  "Dose" = as.numeric(as.vector(sensitivityRaw(pSets[[i]])[
                    exp,
                    ,
                    "Dose"
                  ])),
                  "Viability" = as.numeric(as.vector(sensitivityRaw(pSets[[i]])[
                    exp,
                    ,
                    "Viability"
                  ]))
                ),
                stringsAsFactors = FALSE
              )
              drug.responses <- drug.responses[complete.cases(drug.responses), ]
              doses[[j]] <- drug.responses$Dose
              responses[[j]] <- drug.responses$Viability
              names(doses[[j]]) <- names(responses[[
                j
              ]]) <- seq_len(length(doses[[j]]))
              if (!missing(legends.label)) {
                if (length(legends.label) > 1) {
                  legend.values[[j]] <- paste(
                    unlist(lapply(legends.label, function(x) {
                      sprintf(
                        "%s = %s",
                        x,
                        round(
                          as.numeric(sensitivityProfiles(pSets[[i]])[exp, x]),
                          digits = 2
                        )
                      )
                    })),
                    collapse = ", "
                  )
                } else {
                  legend.values[[j]] <- sprintf(
                    " Exp %s %s = %s",
                    rownames(sensitivityInfo(pSets[[i]]))[exp],
                    legends.label,
                    round(
                      as.numeric(sensitivityProfiles(pSets[[i]])[
                        exp,
                        legends.label
                      ]),
                      digits = 2
                    )
                  )
                }
              } else {
                tt <- unlist(strsplit(
                  rownames(sensitivityInfo(pSets[[i]]))[exp],
                  split = "_"
                ))
                if (tt[1] == treatment_col) {
                  legend.values[[j]] <- tt[2]
                } else {
                  legend.values[[j]] <- rownames(sensitivityInfo(pSets[[i]]))[
                    exp
                  ]
                }
              }
            }
          }
        } else {
          warning(
            sprintf(
              "The cell line and drug combo were not tested together in %s. Skipping.",
              name(pSets[[i]])
            )
          )
          next
        }
      }
    }

    if (!missing(concentrations)) {
      doses2 <- list()
      responses2 <- list()
      legend.values2 <- list()
      j <- 0
      pSetNames2 <- list()
      for (i in seq_len(length(concentrations))) {
        doses2[[i]] <- concentrations[[i]]
        responses2[[i]] <- viabilities[[i]]
        legend_label <- character(0)
        if (length(legends.label) > 0) {
          if (any(grepl("AUC", x = toupper(legends.label)))) {
            auc_label <- sprintf(
              "%s = %s",
              "AUC",
              round(
                computeAUC(
                  concentrations[[i]],
                  viabilities[[i]],
                  conc_as_log = FALSE,
                  viability_as_pct = TRUE
                ) /
                  100,
                digits = 2
              )
            )
            legend_label <- c(legend_label, auc_label)
          }
          if (any(grepl("IC50", x = toupper(legends.label)))) {
            ic50_label <- sprintf(
              "%s = %s",
              "IC50",
              round(
                computeIC50(
                  concentrations[[i]],
                  viabilities[[i]],
                  conc_as_log = FALSE,
                  viability_as_pct = TRUE
                ),
                digits = 2
              )
            )
            legend_label <- c(legend_label, ic50_label)
          }
          legend.values2[[i]] <- paste(legend_label, collapse = ", ")
        } else {
          legend.values2[[i]] <- ""
        }

        pSetNames2[[i]] <- names(concentrations)[[i]]
      }
      doses <- c(doses, doses2)
      responses <- c(responses, responses2)
      legend.values <- c(legend.values, legend.values2)
      pSetNames <- c(pSetNames, pSetNames2)
    }

    if (missing(mycol)) {
      # require(RColorBrewer) || stop("Library RColorBrewer is not available!")
      mycol <- RColorBrewer::brewer.pal(n = 7, name = "Set1")
    }

    dose.range <- c(Inf, -Inf)
    viability.range <- c(0, 10)
    filtered_doses <- vector("list", length(doses))
    filtered_responses <- vector("list", length(responses))
    for (i in seq_len(length(doses))) {
      dose_vec <- doses[[i]]
      resp_vec <- responses[[i]]
      valid_idx <- is.finite(dose_vec) &
        is.finite(resp_vec) &
        !is.na(dose_vec) &
        !is.na(resp_vec) &
        dose_vec > 0
      filtered_doses[[i]] <- dose_vec[valid_idx]
      filtered_responses[[i]] <- resp_vec[valid_idx]
      if (length(filtered_doses[[i]]) > 0) {
        dose.range <- c(
          min(
            dose.range[1],
            min(filtered_doses[[i]], na.rm = TRUE),
            na.rm = TRUE
          ),
          max(
            dose.range[2],
            max(filtered_doses[[i]], na.rm = TRUE),
            na.rm = TRUE
          )
        )
        viability.range <- c(
          0,
          max(
            viability.range[2],
            max(filtered_responses[[i]], na.rm = TRUE),
            na.rm = TRUE
          )
        )
      }
    }
    if (!is.finite(dose.range[1]) || !is.finite(dose.range[2])) {
      warning("No positive finite doses available for plotting.")
      return(invisible(NULL))
    }
    x1 <- 10^10
    x2 <- 0

    if (length(doses) > 1) {
      common.ranges <- .getCommonConcentrationRange(filtered_doses)

      for (i in seq_len(length(doses))) {
        x1 <- min(x1, min(common.ranges[[i]]))
        x2 <- max(x2, max(common.ranges[[i]]))
      }
    }
    if (!missing(xlim)) {
      dose.range <- xlim
    }
    if (!missing(ylim)) {
      viability.range <- ylim
    }
    if (missing(title)) {
      if (!missing(drug) && !missing(cellline)) {
        title <- sprintf("%s:%s", drug, cellline)
      } else {
        title <- "Drug Dose Response Curve"
      }
    }
    plot(
      NA,
      xlab = "Concentration (uM)",
      ylab = "% Viability",
      axes = FALSE,
      main = title,
      log = "x",
      ylim = viability.range,
      xlim = dose.range,
      cex = cex,
      cex.main = cex.main
    )
    magicaxis::magaxis(
      side = seq_len(2),
      frame.plot = TRUE,
      tcl = -.3,
      majorn = c(5, 3),
      minorn = c(5, 2)
    )
    legends <- NULL
    legends.col <- NULL
    if (length(doses) > 1) {
      rect(
        xleft = x1,
        xright = x2,
        ybottom = viability.range[1],
        ytop = viability.range[2],
        col = rgb(240, 240, 240, maxColorValue = 255),
        border = FALSE
      )
    }

    for (i in seq_len(length(doses))) {
      filtered_dose <- filtered_doses[[i]]
      filtered_resp <- filtered_responses[[i]]

      if (length(filtered_dose) < 2 || length(unique(filtered_resp)) < 2) {
        next
      }

      points(filtered_dose, filtered_resp, pch = 20, col = mycol[i], cex = cex)

      ordered_idx <- order(filtered_dose)
      ordered_dose <- filtered_dose[ordered_idx]
      ordered_resp <- filtered_resp[ordered_idx]

      draw_fitted <- function() {
        fit <- try(
          logLogisticRegression(
            conc = filtered_dose,
            viability = filtered_resp
          ),
          silent = TRUE
        )
        if (inherits(fit, "try-error")) {
          return()
        }
        log10_x_vals <- .getSupportVec(log10(filtered_dose))
        pred <- .Hill(
          log10_x_vals,
          pars = c(
            fit$HS,
            fit$E_inf / 100,
            log10(fit$EC50)
          )
        ) *
          100
        if (trunc) {
          pred <- pmin(pmax(pred, 0), 100)
        }
        lines(
          10^log10_x_vals,
          pred,
          lty = 1,
          lwd = lwd,
          col = mycol[i]
        )
      }

      switch(
        plot.type,
        "Actual" = {
          lines(ordered_dose, ordered_resp, lty = 1, lwd = lwd, col = mycol[i])
        },
        "Fitted" = {
          draw_fitted()
        },
        "Both" = {
          lines(ordered_dose, ordered_resp, lty = 1, lwd = lwd, col = mycol[i])
          draw_fitted()
        }
      )
      legends <- c(legends, sprintf("%s%s", pSetNames[[i]], legend.values[[i]]))
      legends.col <- c(legends.col, mycol[i])
    }

    legend(
      legend.loc,
      legend = legends,
      col = legends.col,
      bty = "n",
      cex = cex,
      pch = c(15, 15)
    )
    return(invisible(NULL))
  }
