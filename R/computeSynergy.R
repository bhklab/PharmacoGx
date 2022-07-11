# ==== Loewe Additivity

#' @title Inverse function of Hill equation
#'
#' @description
#' For the dose-response Hill equation of a drug defined by
#' \eqn{E(x) = E_{inf}+\frac{1-E_{inf}}{1+(\frac{x}{EC50})^(\frac{1}{HS})}},
#' that computes the response in viability from a dose in micromole fir a drug,
#' this function is the inverse function of the Hill curve that
#' computes the dose required to produce a given response:
#' \eqn{
#'     f^{-1}(E) = EC50 (
#'     \frac{1-E}{E-E_{inf}} )^{\frac{1}{HS}}
#'     )
#' }
#'
#' @param viability `numeric` is a vector whose entries are the viability values
#'     in the range \[0, 1\].
#' @param EC50 `numeric` is a vector of relative EC50 for drug-response equation.
#' @param HS `numeric` Hill coefficient of the drug-response equation
#'     that represents the sigmoidity of the curve.
#' @param E_inf `numeric` the maximum attanable effect of a drug
#'     when it is administered with a infinitely high concentration.
#'
#' @return `numeric` concentrations in micromoles required to produce
#'     `viability` in the corresponding entries.
#'
#' @examples
#' print("TODO::")
#'
#' @export
effectToDose <- function(viability, EC50, HS, E_inf) {
    ## TODO:: Check input validity
    EC50 * (
        (1 - viability) / (viability - E_inf)
    )^ (1 / HS)
}

#' @title Loewe Additive Combination Index (CI)
#'
#' @description
#' Computes the Loewe additive combination index (CI) from its definition
#' \eqn{
#'     CI = \frac{x_1}{f_1^{-1}(E)} +
#'          \frac{x_2}{f_2^{-1}(E)}
#' }
#'
#' @param viability `numeric` is a vector whose entries are the viability values
#'     in the range \[0, 1\].
#' @param treatment1dose `numeric` a vector of concentrations for treatment 1
#' @param HS_1 `numeric` Hill coefficient of treatment 1
#' @param E_inf_1 `numeric` the maximum attainable effect of treatment 1.
#' @param EC50_1 `numeric` relative EC50 of treatment 1.
#' @param treatment2dose `numeric` a vector of concentrations for treatment 2
#' @param HS_2 `numeric` Hill coefficient of treatment 2
#' @param E_inf_2 `numeric` the maximum attainable effect of treatment 2.
#' @param EC50_2 `numeric` relative EC50 of treatment 2.
#'
#' @return CI under Loewe additive definition
#'
#' @examples
#' print("TODO::")
#'
#' @export
LoeweCI <- function(viability,
                    treatment1dose, HS_1, E_inf_1, EC50_1,
                    treatment2dose, HS_2, E_inf_2, EC50_2) {
    (treatment1dose / effectToDose(
        viability = viability,
        EC50 = EC50_1,
        HS = HS_1,
        E_inf = E_inf_1)) +
    (treatment2dose / effectToDose(
        viability = viability,
        EC50 = EC50_2,
        HS = HS_2,
        E_inf = E_inf_2))
}
## Objective function to mimimise for solving E_Loewe
#' @param viability `numeric` is a vector whose entries are the viability values
#'     in the range [0, 1].
#' @param treatment1dose `numeric` a vector of concentrations for treatment 1
#' @param HS_1 `numeric` Hill coefficient of treatment 1
#' @param E_inf_1 `numeric` the maximum attainable effect of treatment 1.
#' @param EC50_1 `numeric` relative EC50 of treatment 1.
#' @param treatment2dose `numeric` a vector of concentrations for treatment 2
#' @param HS_2 `numeric` Hill coefficient of treatment 2
#' @param E_inf_2 `numeric` the maximum attainable effect of treatment 2.
#' @param EC50_2 `numeric` relative EC50 of treatment 2.
#'
#' @return the distance between computed Loewe CI and 1
#'
#' @noRd
.LoeweLoss <- function(viability,
                       treatment1dose, HS_1, E_inf_1, EC50_1,
                       treatment2dose, HS_2, E_inf_2, EC50_2) {
    abs(
        LoeweCI(viability = viability,
                treatment1dose, HS_1, E_inf_1, EC50_1,
                treatment2dose, HS_2, E_inf_2, EC50_2) - 1
    )
}

#' @title Computes Loewe Null References
#'
#' @description
#' Predict the response of a treatment combination under
#' the Loewe additive null assumption.
#'
#' @param treatment1dose `numeric` a vector of concentrations for treatment 1
#' @param HS_1 `numeric` Hill coefficient of treatment 1
#' @param E_inf_1 `numeric` viability produced by the maximum attainable effect of treatment 1.
#' @param EC50_1 `numeric` relative EC50 of treatment 1.
#' @param treatment2dose `numeric` a vector of concentrations for treatment 2
#' @param HS_2 `numeric` Hill coefficient of treatment 2
#' @param E_inf_2 `numeric` viability produced by the maximum attainable effect of treatment 2.
#' @param EC50_2 `numeric` relative EC50 of treatment 2.
#' @param tol `numeric` Error tolerance for deviations from Loewe assumption. Loewe predictions with error higher than `tol` will be returned as `NA`. Deafult 0.5.
#' @param lower_bound `numeric` Lowest possible value for Loewe expected viability. Default 0.
#' @param upper_bound `numeric` Highest possible value for Loewe expected viability. Default 1.
#'
#' @return `numeric` expected viability under Loewe additive null assumption.
#'
#' @export
#'
#' @examples
#' print("TODO::")
#'
#' @importFrom stats optimise
computeLoewe <- function(treatment1dose, HS_1, E_inf_1, EC50_1,
                         treatment2dose, HS_2, E_inf_2, EC50_2,
                         tol = 0.5, lower_bound = 0, upper_bound = 1) {
    ## Find viability that minimises the distance between Loewe CI and 1
    loewe_guess <- optimise(
        f = .LoeweLoss,
        lower = lower_bound,
        upper = upper_bound,
        treatment1dose = treatment1dose,
        HS_1 = HS_1, E_inf_1 = E_inf_1, EC50_1 = EC50_1,
        treatment2dose = treatment2dose,
        HS_2 = HS_2, E_inf_2 = E_inf_2, EC50_2 = EC50_2
    )

    guess_err <- loewe_guess$objective
    loewe_estimate <- loewe_guess$minimum

    if(guess_err > tol)
        loewe_estimate <- NA_real_

    return(loewe_estimate)
}


# ==== Zero Interaction Potency (ZIP)

#' @title Computes ZIP Null References
#'
#' @description
#' Predict the additive response of a treatment combination under
#' the ZIP null assumption.
#'
#' @param treatment1dose `numeric` a vector of concentrations for treatment 1
#' @param HS_1 `numeric` Hill coefficient of treatment 1
#' @param EC50_1 `numeric` relative EC50 of treatment 1.
#' @param E_inf_1 `numeric` viability produced by the maximum attainable effect of treatment 1.
#'     Default 0 by the original paper.
#' @param treatment2dose `numeric` a vector of concentrations for treatment 2
#' @param HS_2 `numeric` Hill coefficient of treatment 2
#' @param EC50_2 `numeric` relative EC50 of treatment 2.
#' @param E_inf_2 `numeric` viability produced by maximum effect of treatment 2.
#'     Default 0 by the original paper.
#'
#' @return `numeric` expected viability under ZIP null assumption.
#'
#' @export
computeZIP <- function(treatment1dose, HS_1, EC50_1, E_inf_1 = 0,
                       treatment2dose, HS_2, EC50_2, E_inf_2 = 0) {
    y_1 <- .Hill(log10(treatment1dose), c(HS_1, E_inf_1, log10(EC50_1)))
    y_2 <- .Hill(log10(treatment2dose), c(HS_2, E_inf_2, log10(EC50_2)))
    y_zip <- y_1 * y_2
    return(y_zip)
}

#' @title Compute response of adding one treatment to the other
#'
#' @description
#' Response of adding one drug to the other,
#' with the assumption that E_inf = 0, and
#' E_min is dictated by the response of the drug added.
#'
#' @param dose_add `numeric` a fixed concentration of the drug added.
#' @param EC50_add `numeric` relative EC50 of the drug added.
#' @param HS_add `numeric` Hill coefficient of the drug added.
#' @param dose_to `numeric` a vector of concentrations of the drug being added to
#' @param HS_proj `numeric` changed Hill coefficient of the drug being added to; the projected shape parameter.
#' @param EC50_proj `numeric` changed relative EC50 of the drug being added to; the projected potency.
#'
#' @return `numeric` viability after adding one drug to the other.
#'
#' @export
projectedResponse <- function(dose_add,
                              EC50_add,
                              HS_add,
                              dose_to,
                              EC50_proj,
                              HS_proj) {
    E_min <- .Hill(log10(dose_add), c(HS_add, 0, log10(EC50_add)))
    return(
           E_min / (1 + (dose_to/EC50_proj)^(HS_proj))
    )
}

#' L2-loss to optimise for EC50_proj and HS_proj
#'
#' @param par `numeric` EC_proj and HS_proj; the projected/shifted potency and Hill coefficient.
#' @param dose_to `numeric` a vector of concentrations of the drug being added to
#' @param viability `numeric` Observed viability of two treatments; target for fitting curve.
#' @param dose_add `numeric` a vector of concentrations of the drug added.
#' @param HS_add `numeric` Hill coefficient of the drug added.
#' @param EC50_add `numeric` relative EC50 of treatment 1 of the drug added.
#'
#' @return `numeric` L2 loss for fitting a 2-parameter curve. See below.
#'
#' @noRd
.potencyFitLoss <- function(par, dose_to, viability,
                            dose_add, EC50_add, HS_add) {
    ## L2 Loss
    norm(
         projectedResponse(dose_add, EC50_add, HS_add,
                           dose_to,
                           EC50_proj = par[1],
                           HS_proj = par[2]) -
         viability, "2"
    )

}

#' @title Estimate projected potency and shape parameter
#'
#' @description
#' Estimate the projected potency EC50 and the shape parameter HS
#' in the dose-response curve of a drug after adding another drug to it
#' by fitting a 2-parameter dose-response curve.
#' It assumes \eqn{E_min = 1} for the drug being added and
#' \eqn{E_inf = 0} for both drugs.
#'
#' @param dose_to `numeric` a vector of concentrations of the drug being added to
#' @param viability `numeric` Observed viability of two treatments; target for fitting curve.
#' @param dose_add `numeric` a vector of concentrations of the drug added.
#' @param EC50_add `numeric` relative EC50 of the drug added.
#' @param HS_add `numeric` Hill coefficient of the drug added.
#' @param use_L2 `logical` Whether to use L2 loss for fitting curves.
#'     This method produces cruder estimates for projected Hill parameters of
#'     drug combinations difficult to optimise, but also faster to compute results.
#'     Use it only if the default method is not progressing. Default `FALSE`.
#'
#' @return `list` \itemize{
#'      \item{"EC50_proj"}{Projected potency after adding a drug}
#'      \item{"HS_proj"}{Projected Hill coefficient after adding a drug}
#' }
#'
#' @export
#'
#' @importFrom CoreGx .fitCurve
#' @importFrom stats optimise
estimateNewPotency <- function(dose_to, viability,
                               dose_add,
                               EC50_add,
                               HS_add,
                               use_L2 = FALSE) {
    lower_bounds <- c(1e-6, 0)
    upper_bounds <- c(1e+6, 4)
    gritty_guess <- c(
            pmin(pmax(dose_to[which.min(abs(viability - 1/2))],
                      lower_bounds[1]),
                 upper_bounds[1]),
            pmin(pmax(1, lower_bounds[2]), upper_bounds[2])
    )

    if (use_L2) {
        proj_params <- optim(
            fn = .potencyFitLoss,
            par = gritty_guess,
            lower = lower_bounds, upper = upper_bounds,
            method = "L-BFGS-B",
            dose_to = dose_to,
            viability = viability,
            dose_add = dose_add,
            EC50_add = EC50_add,
            HS_add = HS_add
        )$par
    } else {
        ## Default method
        E_min <- .Hill(log10(dose_add), c(HS_add, 0, log10(EC50_add)))
        proj_params <- CoreGx::.fitCurve(
            x = dose_to, y = viability, f = function(dose_to, par) {
                # par[1] = EC50_proj
                # par[2] = HS_proj
                E_min / (1 + (dose_to/par[1])^(par[2]))
            },
            lower_bounds = lower_bounds,
            upper_bounds = upper_bounds,
            density = c(2, 5),
            step = .5 / c(2, 5),
            precision = 1e-4,
            scale = 0.07,
            family = c("normal", "Cauchy"),
            median_n = 1,
            trunc = TRUE,
            verbose = TRUE,
            gritty_guess = gritty_guess
        )
    }

    return(list(
        EC50_proj = proj_params[1],
        HS_proj = proj_params[2]
    ))
}

#' Generic to compute ZIP delta scores from an S4 object
#'
#' @examples
#' print("Generics shouldn't need examples?")
#'
#' @param object `S4` An object to compute delta scores from.
#' @param ... Allow new arguments to this generic.
#'
#' @return Depends on the implemented method.
#'
#' @export
setGeneric(name = "computeZIPdelta",
           def = function(object, ...) standardGeneric("computeZIPdelta"))

#' @title Compute ZIP delta score
#'
#' @description
#' Following the calculation of ZIP delta score as in Appendix A3.
#' See reference for details.
#'
#' @description Compute ZIP delta score as described in the original paper.
#'
#' @param object [TreatmentResponseExperiment]
#'     The `TreatmentResponseExperiment` from which to extract assays
#'     `mono_profile` and `combo_viability` to compute ZIP delta scores.
#' @param use_L2 `logical`
#'     whether to use L2-loss rather than the default optimisation method.
#'     This method produces cruder estimates for projected Hill parameters of
#'     drug combinations difficult to optimise, but also faster to compute results.
#'     Use it only if the default method is not progressing. Default `FALSE`.
#' @param nthread `integer` Number of cores used to perform computation.
#'     Default 1.
#'
#' @return [TreatmentResponseExperiment] with assay `combo_scores` containing `delta_scores`
#'
#' @references
#' Yadav, B., Wennerberg, K., Aittokallio, T., & Tang, J. (2015). Searching for Drug Synergy in Complex Dose–Response Landscapes Using an Interaction Potency Model. Computational and Structural Biotechnology Journal, 13, 504–513. https://doi.org/10.1016/j.csbj.2015.09.001
#'
#' @importFrom CoreGx buildComboProfiles aggregate
#' @import data.table
#' @export
#' @docType methods
setMethod(f = "computeZIPdelta",
          signature = signature(object = "TreatmentResponseExperiment"),
          definition = function(object, use_L2 = FALSE, nthread = 1L) {

    if (!is.logical(use_L2)) {
        stop(.errorMsg("argument `use_L2` must be type of logical"))
    } else if (length(use_L2) != 1) {
        stop(.errorMsg("argument `use_L2` must be of length 1"))
    }

    if (!is.integer(nthread)) {
        stop(.errorMsg("argument `nthread` must be type of integer"))
    } else if (length(nthread) != 1) {
        stop(.errorMsg("argument `nthread` must be of length 1"))
    }

    combo_keys <- c("treatment1id", "treatment2id",
                    "treatment1dose", "treatment2dose", "sampleid")
    combo_profiles <- buildComboProfiles(object, c("HS", "EC50", "E_inf", "viability"))
    setkeyv(combo_profiles, combo_keys)
    combo_profiles |>
        aggregate(
            estimateNewPotency(
                dose_to = treatment1dose,
                viability = viability_combo,
                dose_add = unique(treatment2dose),
                EC50_add = unique(EC50_2),
                HS_add = unique(HS_2),
                use_L2 = use_L2
            ),
            by = c("treatment1id", "treatment2id", "treatment2dose", "sampleid"),
            nthread = nthread,
            enlist = FALSE,
            moreArgs=list(use_L2=use_L2)
        ) -> fit_2_to_1
    combo_profiles |>
        aggregate(
            estimateNewPotency(
                dose_to = treatment2dose,
                viability = viability,
                dose_add = unique(treatment1dose),
                EC50_add = unique(EC50_1),
                HS_add = unique(HS_1),
                use_L2 = use_L2,
            ),
            by = c("treatment1id", "treatment2id", "treatment1dose", "sampleid"),
            nthread = nthread,
            enlist = FALSE,
            moreArgs=list(use_L2=use_L2)
        ) -> fit_1_to_2

    combo_profiles <- combo_profiles[
        fit_1_to_2, ,
        on = c(
            treatment1id = "treatment1id",
            treatment2id = "treatment2id",
            treatment1dose = "treatment1dose",
            sampleid = "sampleid"
        )
    ]

    combo_profiles <- merge.data.table(
        combo_profiles,
        fit_2_to_1,
        by.x = c("treatment1id", "treatment2id", "treatment2dose", "sampleid"),
        by.y = c("treatment1id", "treatment2id", "treatment2dose", "sampleid"),
        suffixes = c("_1_to_2", "_2_to_1")
    )

    combo_profiles |>
        aggregate(
            delta_score = .computeZIPDelta(
                EC50_1_to_2 = EC50_proj_1_to_2,
                EC50_2_to_1 = EC50_proj_2_to_1,
                EC50_1 = EC50_1, EC50_2 = EC50_2,
                HS_1_to_2 = HS_proj_1_to_2,
                HS_2_to_1 = HS_proj_2_to_1,
                HS_1 = HS_1, HS_2 = HS_2,
                E_inf_1 = E_inf_1, E_inf_2 = E_inf_2,
                treatment1dose = treatment1dose,
                treatment2dose = treatment2dose
            ),
            by = combo_keys,
            nthread = nthread
        ) -> delta_scores
    ## Add delta scores to combo_scores in input TRE
    combo_scores <- object$combo_scores
    object$combo_scores <- combo_scores[
        delta_scores, ,
        on = c(treatment1id = "treatment1id",
               treatment2id = "treatment2id",
               treatment1dose = "treatment1dose",
               treatment2dose = "treatment2dose",
               sampleid = "sampleid")]
    ## Can we do call by reference?
    return(object)
})

deltaZIP <- function(treatment1id, treatment2id, treatment1dose, treatment2dose,
        HS_1, HS_2, EC50_1, EC50_2, E_inf_1, E_inf_2) {
    combo_prof <- data.table(
        treatment1id=treatment1id,
        treatment2id=treatment2id,
        treatment1dose=treatment1dose,
        treatment2dose=treatment2dose,
        HS_1=HS_1,
        HS_2=HS_2,
        EC50_1=EC50_1,
        EC50_2=EC50_2,
        E_inf_1=E_inf_1,
        E_inf_2=E_inf_2
    )

    combo_keys <- c("treatment1id", "treatment2id",
                    "treatment1dose", "treatment2dose", "sampleid")
    setkeyv(combo_profiles, combo_keys)

    combo_profiles |>
        aggregate(
            estimateNewPotency(
                dose_to = treatment1dose,
                viability = viability,
                dose_add = unique(treatment2dose),
                EC50_add = unique(EC50_2),
                HS_add = unique(HS_2),
                use_L2 = use_L2
            ),
            by = c("treatment1id", "treatment2id", "treatment2dose", "sampleid"),
            nthread = nthread,
            enlist = FALSE
        ) -> fit_2_to_1
    combo_profiles |>
        aggregate(
            estimateNewPotency(
                dose_to = treatment2dose,
                viability = viability,
                dose_add = unique(treatment1dose),
                EC50_add = unique(EC50_1),
                HS_add = unique(HS_1),
                use_L2 = use_L2
            ),
            by = c("treatment1id", "treatment2id", "treatment1dose", "sampleid"),
            nthread = nthread,
            enlist = FALSE
        ) -> fit_1_to_2

    combo_profiles <- combo_profiles[
        fit_1_to_2, ,
        on = c(
            treatment1id = "treatment1id",
            treatment2id = "treatment2id",
            treatment1dose = "treatment1dose",
            sampleid = "sampleid"
        )
    ]

    combo_profiles <- merge.data.table(
        combo_profiles,
        fit_2_to_1,
        by.x = c("treatment1id", "treatment2id", "treatment2dose", "sampleid"),
        by.y = c("treatment1id", "treatment2id", "treatment2dose", "sampleid"),
        suffixes = c("_1_to_2", "_2_to_1")
    )

    combo_profiles |>
        aggregate(
            delta_score = .computeZIPDelta(
                EC50_1_to_2 = EC50_proj_1_to_2,
                EC50_2_to_1 = EC50_proj_2_to_1,
                EC50_1 = EC50_1, EC50_2 = EC50_2,
                HS_1_to_2 = HS_proj_1_to_2,
                HS_2_to_1 = HS_proj_2_to_1,
                HS_1 = HS_1, HS_2 = HS_2,
                E_inf_1 = E_inf_1, E_inf_2 = E_inf_2,
                treatment1dose = treatment1dose,
                treatment2dose = treatment2dose
            ),
            by = combo_keys,
            nthread = nthread
        ) -> delta_scores

    return(delta_scores$delta_score)
}


#' @title Calculate ZIP delta score for each drug combination
#'
#' @description
#' Following the calculation of ZIP delta score as in Appendix A3.
#' See reference for details.
#'
#' @param EC50_1_to_2 `numeric` projected EC50 of treatment 2 after adding treatment 1.
#' @param EC50_2_to_1 `numeric` projected EC50 of treatment 1 after adding treatment 2.
#' @param EC50_1 `numeric` relative EC50 of treatment 1.
#' @param EC50_2 `numeric` relative EC50 of treatment 2.
#' @param HS_1_to_2 `numeric` projected Hill coefficient of treatment 2
#'     after adding treatment 1.
#' @param HS_2_to_1 `numeric` projected Hill coefficient of treatment 1
#'     after adding treatment 2.
#' @param HS_1 `numeric` Hill coefficient of treatment 1
#' @param HS_2 `numeric` Hill coefficient of treatment 2
#' @param E_inf_1 `numeric` viability produced by the maximum attainable effect of treatment 1.
#' @param E_inf_2 `numeric` viability produced by the maximum attainable effect of treatment 2.
#' @param treatment1dose `numeric` a vector of concentrations for treatment 1
#' @param treatment2dose `numeric` a vector of concentrations for treatment 2
#'
#' @return `numeric` ZIP delta score to quantify synergy.
#' @noRd
#' @references
#' Yadav, B., Wennerberg, K., Aittokallio, T., & Tang, J. (2015). Searching for Drug Synergy in Complex Dose–Response Landscapes Using an Interaction Potency Model. Computational and Structural Biotechnology Journal, 13, 504–513. https://doi.org/10.1016/j.csbj.2015.09.001
#' @export

.computeZIPDelta <- function(EC50_1_to_2, EC50_2_to_1, EC50_1, EC50_2,
                             HS_1_to_2, HS_2_to_1, HS_1, HS_2,
                             E_inf_1, E_inf_2, treatment1dose, treatment2dose) {
    E_min <- 1 ## 3-parameter case; subject to change
    viability_1 <- .Hill(log10(treatment1dose), c(HS_1, E_inf_1, log10(EC50_1)))
    viability_2 <- .Hill(log10(treatment2dose), c(HS_2, E_inf_2, log10(EC50_2)))
    viability_2_to_1 <- (
        viability_2 + E_inf_1*(treatment1dose/EC50_2_to_1)^(HS_2_to_1)
        ) / (
        1 + (treatment1dose/EC50_2_to_1)^(HS_2_to_1)
    )
    viability_1_to_2 <- (
        viability_1 + E_inf_2*(treatment2dose/EC50_1_to_2)^(HS_1_to_2)
        ) / (
        1 + (treatment2dose/EC50_1_to_2)^(HS_1_to_2)
    )
    viability_ZIP <- computeZIP(treatment1dose = treatment1dose,
                                treatment2dose = treatment2dose,
                                HS_1 = HS_1, HS_2 = HS_2,
                                EC50_1 = EC50_1, EC50_2 = EC50_2,
                                E_inf_1 = E_inf_1, E_inf_2 = E_inf_2)
    delta <- (1/2) * (viability_2_to_1 + viability_1_to_2) - viability_ZIP

    return(delta)
}


# ==== Bliss Independence

#' @title Compute Bliss Null References
#'
#' @description
#' Given two `numeric` containing viability of two monotherapy respectively,
#' Compute Bliss null reference values for the expected response
#' of the two treatments combined.
#'
#' @param viability_1 `numeric` monotherapeutic response of treatment 1.
#' @param viability_2 `numeric` monotherapeutic response of treatment 2.
#'
#' @return `numeric` expected response of the two treatments combined
#'     under Bliss null assumption.
#'
#' @export
computeBliss <- function(viability_1, viability_2) {
    bliss_ref <- (viability_1 * viability_2)
    return(bliss_ref)
}

# ==== Highest Single Agent (HSA)

#' @title Compute HSA Null References
#'
#' @description
#' Given two `numeric` containing viability of two monotherapy respectively,
#' Compute highest single-agent (HSA) values as the expected response
#' of the two treatments combined.
#'
#' @param viability_1 `numeric` monotherapeutic response of treatment 1.
#' @param viability_2 `numeric` monotherapeutic response of treatment 2.
#'
#' @return `numeric` expected response of the two treatments combined
#'     using the highest response of the two (lower viability).
#'
#' @export
computeHSA <- function(viability_1, viability_2) {
    HSA_ref <- min(viability_1, viability_2)
    return(HSA_ref)
}
