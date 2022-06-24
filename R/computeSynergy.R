# ==== Loewe Additivity

#' @title Inverse function of Hill equation
#'
#' @description
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
#' @noRd
#'
#' @examples
#' print("TODO::")
#'
#' @export
.effectToDose <- function(viability, EC50, HS, E_inf) {
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
#' @examples
#' print("TODO::")
#'
#' @importFrom stats optimise
computeLoewe <- function(treatment1dose, HS_1, E_inf_1, EC50_1,
                         treatment2dose, HS_2, E_inf_2, EC50_2) {
    ## Find viability that minimises the distance between Loewe CI and 1
    loewe_ref <- optimise(
        f = .LoeweLoss,
        lower = 0,
        upper = 1,
        treatment1dose = treatment1dose,
        HS_1 = HS_1, E_inf_1 = E_inf_1, EC50_1 = EC50_1,
        treatment2dose = treatment2dose,
        HS_2 = HS_2, E_inf_2 = E_inf_2, EC50_2 = EC50_2
    )

    return(loewe_ref$minimum)
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
#' @param treatment2dose `numeric` a vector of concentrations for treatment 2
#' @param HS_2 `numeric` Hill coefficient of treatment 2
#' @param EC50_2 `numeric` relative EC50 of treatment 2.
#' 
#' @return `numeric` expected viability under ZIP null assumption.
#' 
#' @export
computeZIP <- function(treatment1dose,
                       treatment2dose,
                       HS_1, 
                       HS_2, 
                       EC50_1,
                       EC50_2) {
    y_1 <- .Hill(log10(treatment1dose), c(HS_1, 0, log10(EC50_1)))
    y_2 <- .Hill(log10(treatment2dose), c(HS_2, 0, log10(EC50_2)))
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
#' @param dose_add `numeric` a vector of concentrations of the drug added.
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

# Objective function to optimise for EC50_proj and HS_proj
#' @param par `numeric` EC_proj and HS_proj; the projected/shifted potency and Hill coefficient.
#' @param viability `numeric` Observed viability of two treatments; target for fitting curve.
#' @param dose_add `numeric` a vector of concentrations of the drug added.
#' @param HS_add `numeric` Hill coefficient of the drug added.
#' @param EC50_add `numeric` relative EC50 of treatment 1 of the drug added.
#' @param dose_to `numeric` a vector of concentrations of the drug being added to
#' 
#' @return `numeric` L2 loss for fitting a 2-parameter curve. See below.
#' 
#' @noRd
.potencyFitLoss <- function(par, viability,
                            dose_add, EC50_add, HS_add,
                            dose_to) {
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
#' @param viability `numeric` Observed viability of two treatments; target for fitting curve.
#' @param dose_add `numeric` a vector of concentrations of the drug added.
#' @param EC50_add `numeric` relative EC50 of the drug added.
#' @param HS_add `numeric` Hill coefficient of the drug added.
#' @param dose_to `numeric` a vector of concentrations of the drug being added to
#' 
#' @return `list` \itemize{
#'      \item{"EC50"}{Projected potency after adding a drug}
#'      \item{"HS"}{Projected Hill coefficient after adding a drug}
#' }
#' 
#' @export
estimateNewPotency <- function(viability,
                               dose_add,
                               EC50_add,
                               HS_add,
                               dose_to) {
    potency <- optim(
        fn = .potencyFitLoss,
        par = c(1, 0),
        lower = c(1e-6, 0), upper = c(1e+6, 4),
        method = "L-BFGS-B",
        viability = viability,
        dose_add = dose_add,
        EC50_add = EC50_add,
        HS_add = HS_add,
        dose_to = dose_to
    )
    
    return(
        list(
             EC50 = potency$par[1],
             HS = potency$par[2]
        )
    )

}

#' @title Compute ZIP delta score
#'
#' @description Compute ZIP delta score as described in the original paper.
#'
#' @param viability `numeric` Observed viability of two treatments.
#' @param treatment1dose `numeric` a vector of concentrations for treatment 1
#' @param treatment2dose `numeric` a vector of concentrations for treatment 2
#' @param HS_1 `numeric` Hill coefficient of treatment 1
#' @param HS_2 `numeric` Hill coefficient of treatment 2
#' @param EC50_1 `numeric` relative EC50 of treatment 1.
#' @param EC50_2 `numeric` relative EC50 of treatment 2.
#' 
#' @return `numeric` ZIP delta score to quantify synergy.
#' 
#' @export
computeDeltaScore <- function(viability,
                              treatment1dose, treatment2dose,
                              HS_1, HS_2, 
                              EC50_1, EC50_2) {
    potency_1_to_2 <- estimateNewPotency(viability,
                                         dose_add = treatment1dose,
                                         EC50_add = EC50_1,
                                         HS_add = HS_1,
                                         dose_to = treatment2dose)
    potency_2_to_1 <- estimateNewPotency(viability,
                                         dose_add = treatment2dose,
                                         EC50_add = EC50_2,
                                         HS_add = HS_2,
                                         dose_to = treatment1dose)
    y_1_to_2 <- projectedResponse(dose_add = treatment1dose,
                                  EC50_add = EC50_1,
                                  HS_add = HS_1,
                                  dose_to = treatment2dose,
                                  EC50_proj = potency_1_to_2$EC50,
                                  HS_proj = potency_1_to_2$HS)
    y_2_to_1 <- projectedResponse(dose_add = treatment2dose,
                                  EC50_add = EC50_2,
                                  HS_add = HS_2,
                                  dose_to = treatment1dose,
                                  EC50_proj = potency_2_to_1$EC50,
                                  HS_proj = potency_2_to_1$HS)
    y_zip <- computeZIP(treatment1dose = treatment1dose,
                        treatment2dose = treatment2dose,
                        HS_1 = HS_1, 
                        HS_2 = HS_2, 
                        EC50_1 = EC50_1,
                        EC50_2 = EC50_2)
    delta <- (1/2) * (y_2_to_1 + y_1_to_2) - y_zip

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

