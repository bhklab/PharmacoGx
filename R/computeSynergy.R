# ==== Loewe Additivity

#' @title Inverse function of Hill equation
#'
#' For the dose-response Hill equation of a drug defined by
#' \eqn{E(x) = E_{inf}+\frac{1-E_{inf}}{1+(\frac{x}{EC50})^(\frac{1}{HS})}},
#' that computes the response in viability from a dose in micromole fir a drug,
#' this function is the inverse function of the Hill curve that
#' computes the dose required from a response:
#' \eqn{
#'     f^{-1}(E) = EC50 (
#'     \frac{1-E}{E-E_{inf}} )^{\frac{1}{HS}}
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
#' @export
.effectToDose <- function(viability, EC50, HS, E_inf) {
    ## TODO:: Check input validity
    EC50 * (
        (1 - viability) / (viability - E_inf)
    )^ (1 / HS)
}

#' @title Loewe Additive Combination Index (CI)
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
#' @export
.loewe <- function(viability,
                   treatment1dose, HS_1, E_inf_1, EC50_1,
                   treatment2dose, HS_2, E_inf_2, EC50_2) {
    (treatment1dose / .effectToDose(
        viability = viability,
        EC50 = EC50_1,
        HS = HS_1,
        E_inf = E_inf_1)) +
    (treatment2dose / .effectToDose(
        viability = viability,
        EC50 = EC50_2,
        HS = HS_2,
        E_inf = E_inf_2))
}

## Objective function to minimise for solving E_Loewe
.loeweLoss <- function(viability,
                       treatment1dose, HS_1, E_inf_1, EC50_1,
                       treatment2dose, HS_2, E_inf_2, EC50_2) {
    abs(
        .loewe(viability = viability,
               treatment1dose, HS_1, E_inf_1, EC50_1,
               treatment2dose, HS_2, E_inf_2, EC50_2) - 1
    )
}

#' @title Computes Loewe Null References
#'
#' Predict the response of a treatment combination under
#' the Loewe additive null assumption.
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
#' @return `numeric` expected viability under Loewe additive null assumption.
#'
#' @export
#'
#' @importFrom stats optimise
.computeLoewe <- function(viability,
                          treatment1dose, treatment2dose,
                          HS_1, HS_2,
                          E_inf_1, E_inf_2,
                          EC50_1, EC50_2) {
    loewe_ref <- optimise(
        f = .loeweLoss,
        lower = 0,
        upper = 1,
        treatment1dose = treatment1dose,
        HS_1 = HS_1, E_inf_1 = E_inf_1, EC50_1 = EC50_1,
        treatment2dose = treatment2dose,
        HS_2 = HS_2, E_inf_2 = E_inf_2, EC50_2 = EC50_2
    )

    return(loewe_ref$minimum)
}


# ==== Bliss Independence

##' Compute Bliss Independence score from `numeric`
##' Return the Bliss independence score given the effect of a drug combination.
##'
##' @param y1 `numerical`
##' @param y2 `numerical`
##' @param y  `numerical`
##'
##' @return A `list`
##'
##' @export
##' @noRd
#.computeBliss <- function(y1, y2, y, score = "diff") {
#    CI_metrics <- list(
#        ## TODO:: Add Hill slop fit for Bliss estimates
#        "diff" = quote(bliss_ref - y)
#    )
#    ## Since we are using viability instead of cell death
#    bliss_ref <- (y1 * y2)
#    return(
#        list(
#            "Bliss" = bliss_ref,
#            "Bliss_score" = eval(CI_metrics[[CI]])
#        )
#    )
#}
#
## .computeHSA <- function(y1, y2, y) {
##     HSA_ref <- min(y1, y2) ## should be min for viability
##     return(
##         list(
##             "HSA" = HSA_ref,
##             "HSA_score" = HSA_ref - y
##         )
##     )
#
## }
#
#