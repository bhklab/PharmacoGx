# ==== Loewe Additivity

#' @title Inverse function of Hill equation
#'
#' @description
#' For the dose-response Hill equation of a drug defined by
#' \eqn{E(x) = E_{inf}+\frac{1-E_{inf}}{1+(\frac{x}{EC50})^(\frac{1}{HS})}},
#' that computes the response in viability from a dose in micromole for a drug,
#' this function is the inverse function of the Hill curve that
#' computes the dose required to produce a given response:
#' \eqn{
#'     f^{-1}(E) = EC50 (
#'     \frac{1-E}{E-E_{inf}} )^{\frac{1}{HS}}
#'     )
#' }
#'
#' @param viability `numeric` is a vector whose entries are the viability values
#'     in the range \[0, 1\] if `is_pct` is `FALSE` or \[0, 100\] if it is
#'     `TRUE`.
#' @param EC50 `numeric` is a vector of relative EC50 for drug-response equation.
#' @param HS `numeric` Hill coefficient of the drug-response equation
#'     that represents the sigmoidity of the curve.
#' @param E_inf `numeric` the maximum attanable effect of a drug
#'     when it is administered with a infinitely high concentration.
#' @param is_pct `logical` whether both the input viabiliy and `E_inf` are given
#'     in percentage (\[0, 100\]) rather than decimal (\[0, 1\]). Default FALSE.
#'
#' @return `numeric` concentrations in micromoles required to produce
#'     `viability` in the corresponding entries.
#'
#' @examples
#' dose <- effectToDose(viability = 80,
#'                      EC50 = 42,
#'                      HS = 1,
#'                      E_inf = 10,
#'                      is_pct = TRUE)
#'
#' @importFrom checkmate assertLogical
#' @export
effectToDose <- function(viability, EC50, HS, E_inf, is_pct = FALSE) {
    assertLogical(is_pct, len = 1)
    if (is_pct) {
        viability <- viability / 100
        E_inf <- E_inf / 100
    }
    EC50 * ((1 - viability) / (viability - E_inf))^(1 / HS)
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
#' @param is_pct `logical` whether both the input viabiliy and E_inf are given
#'     in percentage (\[0, 100\]) rather than decimal (\[0, 1\]). Default FALSE.
#'
#' @return CI under Loewe additive definition
#'
#' @examples
#' \dontrun{
#' tre |>
#'     endoaggregate(
#'         assay="combo_viability",
#'         Loewe = PharmacoGx::computeLoewe(
#'             treatment1dose = treatment1dose,
#'             treatment2dose = treatment2dose,
#'             HS_1 = HS_1,
#'             HS_2 = HS_2,
#'             E_inf_1 = E_inf_1,
#'             E_inf_2 = E_inf_2,
#'             EC50_1 = EC50_1,
#'             EC50_2 = EC50_2
#'         ),
#'         by = assayKeys(tre, "combo_viability")
#'     ) -> tre
#' }
#'
#' @export
loeweCI <- function(viability,
                    treatment1dose, HS_1, E_inf_1, EC50_1,
                    treatment2dose, HS_2, E_inf_2, EC50_2,
                    is_pct = FALSE) {
    (treatment1dose / effectToDose(
        viability = viability,
        EC50 = EC50_1,
        HS = HS_1,
        E_inf = E_inf_1,
        is_pct = is_pct)) +
    (treatment2dose / effectToDose(
        viability = viability,
        EC50 = EC50_2,
        HS = HS_2,
        E_inf = E_inf_2,
        is_pct = is_pct))
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
.loeweLoss <- function(viability,
                       treatment1dose, HS_1, E_inf_1, EC50_1,
                       treatment2dose, HS_2, E_inf_2, EC50_2) {
    abs(
        loeweCI(viability = viability,
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
#' @param tol `numeric` Error tolerance for deviations from Loewe assumption. Loewe predictions with error higher than `tol` will be returned as `NA`. Deafult 0.1.
#' @param lower_bound `numeric` Lowest possible value for Loewe expected viability. Default 0.
#' @param upper_bound `numeric` Highest possible value for Loewe expected viability. Default 1.
#' @param verbose `logical` whether to display warning messages. Default `FALSE`.
#'
#' @return `numeric` expected viability under Loewe additive null assumption.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' tre |>
#'     endoaggregate(
#'         assay="combo_viability",
#'         Loewe = computeLoewe(
#'             treatment1dose=treatment1dose,
#'             treatment2dose=treatment2dose,
#'             HS_1=HS_1,
#'             HS_2=HS_2,
#'             E_inf_1=E_inf_1,
#'             E_inf_2=E_inf_2,
#'             EC50_1=EC50_1,
#'             EC50_2=EC50_2
#'         ),
#'         by = assayKeys(tre, "combo_viability")
#'     ) -> tre
#' }
#'
#' @importFrom stats optimise
#' @importFrom checkmate assertNumeric assertLogical
computeLoewe <- function(treatment1dose, HS_1, E_inf_1, EC50_1,
                         treatment2dose, HS_2, E_inf_2, EC50_2,
                         tol = 0.1, lower_bound = 0, upper_bound = 1,
                         verbose = FALSE) {

    len <- length(treatment1dose)
    assertNumeric(treatment1dose, len = len)
    assertNumeric(treatment2dose, len = len)
    assertNumeric(HS_1, len = len)
    assertNumeric(HS_2, len = len)
    assertNumeric(E_inf_1, len = len)
    assertNumeric(E_inf_2, len = len)
    assertNumeric(EC50_1, len = len)
    assertNumeric(EC50_2, len = len)
    assertNumeric(tol, len = 1)
    assertNumeric(lower_bound, len = 1)
    assertNumeric(upper_bound, len = 1)
    assertLogical(verbose, len = 1)

    ## Find viability that minimises the distance between Loewe CI and 1
    if (verbose) {
        loewe_guess <- optimise(
            f = .loeweLoss,
            lower = lower_bound,
            upper = upper_bound,
            treatment1dose = treatment1dose,
            HS_1 = HS_1, E_inf_1 = E_inf_1, EC50_1 = EC50_1,
            treatment2dose = treatment2dose,
            HS_2 = HS_2, E_inf_2 = E_inf_2, EC50_2 = EC50_2
        )
    } else {
        suppressWarnings({
            loewe_guess <- optimise(
                f = .loeweLoss,
                lower = lower_bound,
                upper = upper_bound,
                treatment1dose = treatment1dose,
                HS_1 = HS_1, E_inf_1 = E_inf_1, EC50_1 = EC50_1,
                treatment2dose = treatment2dose,
                HS_2 = HS_2, E_inf_2 = E_inf_2, EC50_2 = EC50_2
            )
        })
    }

    guess_err <- loewe_guess$objective
    loewe_estimate <- loewe_guess$minimum

    if (is.nan(guess_err) | guess_err > tol)
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
#' @examples
#' (zip <- computeZIP(
#'   treatment1dose = c(0.1, 0.01, 0.001),
#'   treatment2dose = c(1, 0.1, 0.01),
#'   HS_1 = rep(1, 3), HS_2 = rep(1.2, 3),
#'   EC50_1 = rep(0.01, 3), EC50_2 = rep(0.1, 3),
#'   E_inf_1 = rep(0, 3), E_inf_2 = rep(0.1, 3)
#' ))
#'
#' @importFrom checkmate assertNumeric
#'
#' @export
computeZIP <- function(treatment1dose, HS_1, EC50_1, E_inf_1,
                       treatment2dose, HS_2, EC50_2, E_inf_2) {
    len <- length(treatment1dose)
    assertNumeric(treatment1dose, len = len)
    assertNumeric(treatment2dose, len = len)
    assertNumeric(HS_1, len = len)
    assertNumeric(HS_2, len = len)
    assertNumeric(E_inf_1, len = len)
    assertNumeric(E_inf_2, len = len)
    assertNumeric(EC50_1, len = len)
    assertNumeric(EC50_2, len = len)

    y_1 <- .Hill(log10(treatment1dose), c(HS_1, E_inf_1, log10(EC50_1)))
    y_2 <- .Hill(log10(treatment2dose), c(HS_2, E_inf_2, log10(EC50_2)))
    y_zip <- y_1 * y_2
    return(y_zip)
}

#' @title 4-Parameter Hill Equation for Stimuli-Response Curves
#'
#' @description
#' Sigmoidal function which fits well to many stimuli-response associations
#' observed in biology and pharmacology. In the context of PharmacoGx we
#' are using it to model treatment-response assocations in cancer cell lines.
#'
#' @param dose `numeric()` A vector of `log10(dose)` values (or equivalent for
#' the stimuli being modelleled).
#' @param HS `numeric(1)` Hill coefficient (n) which defines the slope of the
#' dose-response curve at the mid-point. This parameter describes the degree
#' of sigmoidicity of the Hill curve. HS = 1 corresponds to the rectangular
#' hyperbola in dose-response space.
#' @param EC50 `numeric(1)` The dose required to produce 50% of the
#' theoretically maximal response in the system, `E_inf`. Should be in the same
#' units as `dose`!
#' @param E_inf `numeric(1)` Theoretical maximal response (minimal viability)
#' in the system as a proportion in the range \\[0, 1\\]. Note that since we are
#' predicting viability (percent of cells alive after treatment) instead of
#' response, this value should be low (i.e., more cell killing).
#' @param E_ninf `numeric(1)` Theoretical minimum response (basal response).
#' Defaults to 1, which should be the case for most viability experiments since
#' we expect no cell killing to occur prior to applying a treatment.
#'
#' @return `numeric()` Vector of predicted viabilities for the Hill curve defined
#' by `EC50`, `E_inf`, `E_ninf` and `HS` for each supplied value of `dose`.
#'
#' @references
#' Gesztelyi, R., Zsuga, J., Kemeny-Beke, A., Varga, B., Juhasz, B., &
#' Tosaki, A. (2012). The Hill equation and the origin of quantitative
#' pharmacology. Archive for History of Exact Sciences, 66(4), 427–438.
#' https://doi.org/10.1007/s00407-012-0098-5
#'
#' Motulsky, H., & Christopoulos, A. (2004). Fitting models to biological data
#' using linear and nonlinear regression: A practical guide to curve fitting.
#' Oxford University Press. See Chapter 41.
#'
#' @author
#' Feifei Li
#' Petr Smirnov
#' Christopher Eeles
#'
#' @examples
#' (viability <- hillCurve(
#'   dose=c(0.1, 0.01, 0.001),
#'   HS=1.1,
#'   EC50=0.01,
#'   E_ninf=1,
#'   E_inf=0
#' ))
#'
#' @export
hillCurve <- function(dose, HS, EC50, E_inf, E_ninf) {
    E_inf + (( E_ninf - E_inf ) / ( 1 + ( 10^dose / 10^EC50 )^(HS) ))
}

## TODO:: If it works well for fitting 2-way Hill curves, move it to CoreGx

#' @title Compute Logarithm of Hyperbolic Cosine function
#'
#' @description
#' A numerical stable version of `log(cosh(x))`
#' without floating overflow or underflow.
#' Originally implemented in `limma` by Gordon K Smyth.
#'
#' @param x `numeric` vector or matrix.
#'
#' @return `numeric` a vector or matrix with the same dimension as `x`.
#'
#' @references
#' Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK (2015). “limma powers differential expression analyses for RNA-sequencing and microarray studies.” Nucleic Acids Research, 43(7), e47. doi: 10.1093/nar/gkv007.
#'
#' @noRd
#' @export
.logcosh <- function(x) {
    y <- abs(x) - log(2)
    i <- abs(x) < 1e-4
    y[i] <- 0.5*x[i]^2
    i <- !i & (abs(x) < 17)
    y[i] <- log(cosh(x[i]))
    y
}


#' @title Log-cosh loss for fitting projected Hill curves
#'
#' @description
#' Compute the log hyperbolic cosine (log-cosh) loss,
#' which behaves as L2 at small values and as L1 at large values.
#'
#' @param par `numeric` a vector of parameters to optimise in the following order:
#'     `c(HS_proj, E_inf_proj, EC50_proj)`
#' @param dose_to `numeric` a vector of concentrations of the drug being added to
#' @param viability `numeric` Observed viability of two treatments; target for fitting curve.
#' @param E_min_proj `numeric` Projected `E_min` given by
#'     the viability of the added treatment at a fixed dose.
#'
#' @return `numeric` Log-Cosh loss for fitting a 3-parameter Hill curve. See below.
#'
#' @noRd
.fitProjParamsLoss <- function(par, dose_to, viability, E_min_proj) {
    sum(
        .logcosh(
             hillCurve(
                dose = dose_to,
                E_ninf = E_min_proj,
                HS = par[1],
                E_inf = par[2],
                EC50 = par[3]
            ) - viability
        )
    )
}

#' @title Log-cosh loss for fitting projected Hill curves
#'
#' @description
#' Compute the log hyperbolic cosine (log-cosh) loss,
#' which behaves as L2 at small values and as L1 at large values.
#'
#' @param par `numeric` a vector of parameters to optimise
#' @param x `numeric` a vector of input values to the model
#' @param y `numeric` a vector of target values
#' @param fn `numeric` model to fit
#' @param ... `pairlist` Fall through arguments to `fn`.
#'
#' @return `numeric` scalar Log-Cosh loss for fitting a curve.
#'
#' @keywords interal
#' @noRd
.logcoshLoss <- function(par, x, y, fn, ...) {
    sum(.logcosh(fn(par = par, x) - y))
}

#' @title Estimate the projected Hill coefficient, efficacy, and potency
#'
#' @description
#' Estimate the projected shape parameter HS, efficacy `E_inf` and potency `EC50`
#' in the new dose-response curve of a drug after adding another drug to it
#' by fitting a 2-parameter dose-response curve.
#'
#' @param dose_to `numeric` a vector of concentrations of the drug being added to
#' @param combo_viability `numeric` observed viability of two treatments; target for fitting curve.
#' @param dose_add `numeric` a vector of concentrations of the drug added.
#' @param EC50_add `numeric` relative EC50 of the drug added.
#' @param HS_add `numeric` Hill coefficient of the drug added.
#' @param E_inf_add `numeric` Efficacy of the drug added.
#' @param residual `character` Method used to minimise residual in fitting curves.
#'     3 methods available: `logcosh`, `normal`, `Cauchy`.
#'     The default method is `logcosh`.
#'     It minimises the logarithmic hyperbolic cosine loss of the residuals
#'     and provides the fastest estimation among the three methods,
#'     with fitting quality in between `normal` and `Cauchy`;
#'     recommanded when fitting large-scale datasets.
#'     The other two methods minimise residuals by
#'     considering the truncated probability distribution (as in their names) for the residual.
#'     `Cauchy` provides the best fitting quality but also takes the longest to run.
#' @param show_Rsqr `logical` whether to show goodness-of-fit value in the result.
#' @param conc_as_log `logical` indicates whether input concentrations are in log10 scale.
#' @param loss_args `list` Additional argument to the `loss` function.
#'   These get passed to losss via `do.call` analagously to using `...`.
#' @param optim_only `logical(1)` Should the fall back methods when optim fails
#'
#' @references
#' Motulsky, H., & Christopoulos, A. (2004). Fitting dose-response curves. In Fitting models to biological data using linear and nonlinear regression: A practical guide to curve fitting. Oxford University Press.
#'
#' @return `list`
#'      * `HS_proj`: Projected Hill coefficient after adding a drug
#'      * `E_inf_proj`: Projected efficacy after adding a drug
#'      * `EC50_proj`: Projected potency after adding a drug
#'      * `E_ninf_proj`: Projected baseline viability by the added drug
#'      * `Rsqr`: if `show_Rsqr` is `TRUE`, it will include the R squared value indicating the quality of the fit in the result.
#'
#' @importFrom CoreGx .fitCurve2 .reformatData
#' @importFrom checkmate assertNumeric assertLogical
#'
#' @export
estimateProjParams <- function(dose_to, combo_viability, dose_add, EC50_add, HS_add,
    E_inf_add = 0,
    residual = c("logcosh", "normal", "Cauchy"),
    show_Rsqr = TRUE,
    conc_as_log = FALSE,
    optim_only = FALSE,
    loss_args = list()
) {

    len_to <- length(dose_to)
    assertNumeric(dose_to, len = len_to)
    assertNumeric(combo_viability, len = len_to)
    assertNumeric(dose_add, len = 1)
    assertNumeric(EC50_add, len = 1)
    assertNumeric(HS_add, len = 1)
    assertNumeric(E_inf_add, len = 1)
    assertLogical(show_Rsqr, len = 1)
    assertLogical(conc_as_log, len = 1)
    residual <- match.arg(residual)

    ## viability of the drug being added as the minimum baseline response
    if (conc_as_log) {
        E_ninf_proj <- .Hill(dose_add, c(HS_add, E_inf_add, log10(EC50_add)))
    } else {
        E_ninf_proj <- .Hill(log10(dose_add), c(HS_add, E_inf_add, log10(EC50_add)))
    }
    formatted_data <- .reformatData(
        x = dose_to,
        y = combo_viability,
        x_to_log = !conc_as_log,
        y_to_frac = FALSE, ## subject to change
        y_to_log = FALSE,
        trunc = FALSE
    )
    log_conc <- formatted_data[["x"]]
    combo_viability <- formatted_data[["y"]]

    residual_fns <- list(
        "normal" = CoreGx:::.normal_loss,
        "Cauchy" = CoreGx:::.cauchy_loss,
        "logcosh" = .logcoshLoss
    )
    ## c(HS, EC50, E_inf, E_ninf)
    lower_bounds <- c(0, -6, 0)
    upper_bounds <- c(4, 6, 1)
    density <- c(2, 5, 10)
    step <- 0.5 / density
    gritty_guess <- c(
        pmin(pmax(1, lower_bounds[1]), upper_bounds[1]),
        pmin(
            pmax(
                log_conc[which.min(abs(combo_viability - 1/2))],
                lower_bounds[2]
            ),
            upper_bounds[2]
        ),
        pmin(pmax(min(combo_viability), lower_bounds[3]), upper_bounds[3])
    )

    ## If we have zero or less degrees of freedom, fix the HS parameter to 1
    ## This is as per recommendations in Motulsky & Christopoulos (2004)
    insuff_df <- len_to <= 3
    fit_curve_args <- list(
            par = if (insuff_df) gritty_guess[-1] else gritty_guess,
            x = log_conc,
            y = combo_viability,
            fn = function(x, HS, EC50, E_inf, E_ninf) {
                hillCurve(dose=x, HS, EC50, E_inf, E_ninf)
            },
            loss = residual_fns[[residual]],
            lower = if (insuff_df) lower_bounds[-1] else lower_bounds,
            upper = if (insuff_df) upper_bounds[-1] else upper_bounds,
            density = if(insuff_df) density[-1] else density,
            step = if (insuff_df) step[-1] else step,
            optim_only = optim_only,
            loss_args = loss_args,
            E_ninf = E_ninf_proj
    )
    if (insuff_df)
        fit_curve_args <- c(fit_curve_args, HS = 1)

    proj_params <- do.call(.fitCurve2, fit_curve_args)
    if (insuff_df)
        proj_params <- c(1, proj_params)

    proj_params[2] <- 10^proj_params[2]

    if (show_Rsqr) {
        Rsqr <- attr(proj_params, "Rsquare")
        return(list(
            HS_proj = proj_params[1],
            EC50_proj = proj_params[2],
            E_inf_proj = proj_params[3],
            E_ninf_proj = E_ninf_proj,
            Rsqr = Rsqr
        ))
    } else {
        return(list(
            HS_proj = proj_params[1],
            E_inf_proj = proj_params[2],
            EC50_proj = proj_params[3],
            E_ninf_proj = E_ninf_proj
        ))
    }
}

#' @title Two-way fitting for projected dose-response curve.
#'
#' @description
#' Fit projected dose-response curves with `E_min` as the viability
#' of the treatment being added to the other treament at a fixed dose.
#'
#' @examples
#' \dontrun{
#' combo_profiles <- CoreGx::buildComboProfiles(tre, c("HS", "EC50", "E_inf", "viability"))
#' combo_twowayFit <- fitTwowayZIP(combo_profiles)
#' }
#'
#' @param combo_profiles [data.table] contains three parameters of dose-response curves
#'     for each single agent in a drug comnbination,
#'     and the observed viability of two treatments combined.
#'
#' @param residual `character` Method used to minimise residual in fitting curves.
#'     3 methods available: `c("logcosh", "normal", "Cauchy")`.
#'     The default method is `logcosh`.
#'     It minimises the logarithmic hyperbolic cosine loss of the residuals
#'     and provides the fastest estimation among the three methods,
#'     with fitting quality in between `normal` and `Cauchy`;
#'     recommanded when fitting large-scale datasets.
#'     The other two methods minimise residuals by
#'     considering the truncated probability distribution (as in their names) for the residual.
#'     `Cauchy` provides the best fitting quality but also takes the longest to run.
#' @param show_Rsqr `logical` whether to show goodness-of-fit value in the result.
#' @param nthread `integer` Number of cores used to perform computation. Default 1.
#' @param loss_args `list` Additional argument to the `loss` function.
#'   These get passed to losss via `do.call` analagously to using `...`.
#' @param optim_only `logical(1)` Should the fall back methods when optim fails
#'
#' @return [data.table] contains parameters of projected dose-response curves
#'    for adding one treatment to the other.
#'
#' @references
#' Yadav, B., Wennerberg, K., Aittokallio, T., & Tang, J. (2015). Searching for Drug Synergy in Complex Dose–Response Landscapes Using an Interaction Potency Model. Computational and Structural Biotechnology Journal, 13, 504–513. https://doi.org/10.1016/j.csbj.2015.09.001
#'
#' @importFrom CoreGx aggregate
#' @importFrom checkmate assertLogical assertInt assertDataTable
#' @import data.table
#' @export
fitTwowayZIP <- function(
    combo_profiles,
    residual = "logcosh",
    show_Rsqr = TRUE,
    nthread = 1L,
    optim_only = TRUE,
    loss_args = list()
) {

    assertDataTable(combo_profiles, min.rows = 1)
    assertLogical(show_Rsqr, len = 1)
    assertInt(nthread, lower = 1L)
    required_cols <- c(
        "treatment1id", "treatment2id", "treatment1dose", "treatment2dose",
        "sampleid", "combo_viability",
        "HS_1", "HS_2", "E_inf_1", "E_inf_2", "EC50_1", "EC50_2"
    )

    has_cols <- required_cols %in% colnames(combo_profiles)
    if (!all(has_cols))
        stop("Missing required columns of parameters: ",
             paste(required_cols[!has_cols], sep = ", "),
             call. = FALSE)

    combo_profiles |>
        aggregate(
            estimateProjParams(
                dose_to = treatment1dose,
                combo_viability = combo_viability,
                dose_add = unique(treatment2dose),
                EC50_add = unique(EC50_2),
                HS_add = unique(HS_2),
                E_inf_add = unique(E_inf_2),
                residual = residual,
                show_Rsqr = show_Rsqr,
                optim_only = optim_only,
                loss_args = loss_args
            ),
            moreArgs = list(
                residual = residual,
                show_Rsqr = show_Rsqr,
                optim_only = optim_only,
                loss_args = loss_args
            ),
            by = c("treatment1id", "treatment2id", "treatment2dose", "sampleid"),
            nthread = nthread,
            enlist = FALSE
        ) -> fit_2_to_1
    combo_profiles |>
        aggregate(
            estimateProjParams(
                dose_to = treatment2dose,
                combo_viability = combo_viability,
                dose_add = unique(treatment1dose),
                EC50_add = unique(EC50_1),
                HS_add = unique(HS_1),
                E_inf_add = unique(E_inf_1),
                residual = residual,
                show_Rsqr = show_Rsqr,
                optim_only = optim_only,
                loss_args = loss_args
            ),
            moreArgs = list(
                residual = residual,
                show_Rsqr = show_Rsqr,
                optim_only = optim_only,
                loss_args = loss_args
            ),
            by = c("treatment1id", "treatment2id", "treatment1dose", "sampleid"),
            nthread = nthread,
            enlist = FALSE
        ) -> fit_1_to_2

    combo_twowayFit <- combo_profiles[
        fit_1_to_2, ,
        on = c(
            treatment1id = "treatment1id",
            treatment2id = "treatment2id",
            treatment1dose = "treatment1dose",
            sampleid = "sampleid"
        )
    ]

    combo_twowayFit <- merge.data.table(
        combo_twowayFit,
        fit_2_to_1,
        by.x = c("treatment1id", "treatment2id", "treatment2dose", "sampleid"),
        by.y = c("treatment1id", "treatment2id", "treatment2dose", "sampleid"),
        suffixes = c("_1_to_2", "_2_to_1")
    )

    return(combo_twowayFit)
}

## == Plot the result of two-way fittings for a drug combination experiment ===

#' @title Plot projected Hill curves
#'
#' @description
#' Plot the two-way projected Hill curves of adding one drug to the other.
#'
#' @param combo_twowayFit `data.table`
#'      containing two-way fitted parameters for multiple drug combination experiments.
#'
#' @param treatment1 `character`
#'      the `treatment1id` to select in `combo_twowayFit` for a drug combination.
#'
#' @param treatment2
#'      the `treatment2id` to select in `combo_twowayFit` for a drug combination.
#'
#' @param cellline
#'      the `sampleid` to select in `combo_twowayFit` for a drug combination experiment.
#'
#' @param add_treatment
#'      The added treatment in projected Hill curves, either integer 1 or 2.
#'      1 means adding treatment 1 to treatment 2.
#'
#' @return produce a plot with projected Hill curves of adding treatment [add_treatment]
#'
#' @importFrom graphics plot curve points legend
#' @importFrom grDevices palette rainbow
#' @export
#' @noRd
#' @examples
#' \dontrun{
#' combo_profiles <- CoreGx::buildComboProfiles(tre, c("HS", "EC50", "E_inf", "viability"))
#' combo_twowayFit <- fitTwowayZIP(combo_profiles)
#' .plotProjHill(combo_twowayFit,
#'               treatment1 = "Methotrexate",
#'               treatment2 = "Zolendronic Acid",
#'               cellline = "UO-31",
#'               add_treatment = 1)
#' }
.plotProjHill <- function(combo_twowayFit, treatment1, treatment2,
                          cellline, add_treatment = 1, title = NULL) {

    required_cols <- c("treatment1id", "treatment1dose", "treatment2id", "treatment2dose",
                       "sampleid", "combo_viability", "HS_1", "E_inf_1", "EC50_1", "HS_2", "E_inf_2",
                       "EC50_2", "HS_proj_1_to_2", "E_inf_proj_1_to_2", "EC50_proj_1_to_2",
                       "E_ninf_proj_1_to_2", "HS_proj_2_to_1", "E_inf_proj_2_to_1",
                       "EC50_proj_2_to_1", "E_ninf_proj_2_to_1")
    has_cols <- (required_cols %in% colnames(combo_twowayFit))
    if (!all(has_cols))
        stop("Missing required columns for plotting: ",
             paste(required_cols[!has_cols]))

    select_combo <- combo_twowayFit[treatment1id == treatment1 &
                                    treatment2id == treatment2 &
                                    sampleid == cellline]
    if (dim(select_combo)[1] <= 0)
        stop(paste("No such drug combination with treatment1id:", treatment1,
                   "and treatment2id:", treatment2, "and sampleid:", cellline))

    if (length(add_treatment) > 1)
        stop("Argument `add_treatment` must be of length 1.")

    if (!(add_treatment %in% c(1, 2)))
        stop("Argument `add_treatment` must be either 1 or 2.")

    ## Use variable name as title if not provided
    if (is.null(title))
        title <- deparse(substitute(combo_twowayFit))

    ## Colours for each curve of a fixed concentration of the drug added

    has_Rsqr <- c("Rsqr_1_to_2", "Rsqr_2_to_1") %in% colnames(combo_twowayFit)
    if (add_treatment == 2) {
        ## unique treatment 2 concentrations
        unique_t2_dose <- unique(select_combo[, treatment2dose])
        cols <- palette(rainbow(length(unique_t2_dose)))
        if (has_Rsqr[2])
            Rsqr_2_to_1 <- vector(mode = "numeric", length = length(unique_t2_dose))
        ## Initialise an empty background canvas
        plot(
            NULL, xlim = c(-10, 10), ylim = c(0, 2),
            ylab = paste("Response of adding", treatment2, "to", treatment1),
            xlab = paste0("log10([", treatment1,"])"),
            main = title
        )
        for (i in seq_along(unique_t2_dose)) {
            dose_add <- unique_t2_dose[i]
            EC50_proj <- unique(select_combo[treatment2dose == dose_add, EC50_proj_2_to_1])
            HS_proj <- unique(select_combo[treatment2dose == dose_add, HS_proj_2_to_1])
            EC50_add <- unique(select_combo[treatment2dose == dose_add, EC50_2])
            E_inf_add <- unique(select_combo[treatment2dose == dose_add, E_inf_2])
            E_inf_proj <- unique(select_combo[treatment2dose == dose_add, E_inf_proj_2_to_1])
            HS_add <- unique(select_combo[treatment2dose == dose_add, HS_2])
            dose_to <- select_combo[treatment2dose == dose_add, treatment1dose]
            E_ninf_proj <- PharmacoGx:::.Hill(log10(dose_add), c(HS_add, E_inf_add, log10(EC50_add)))
            if (has_Rsqr[2])
                Rsqr_2_to_1[i] <- unique(select_combo[treatment2dose == dose_add, Rsqr_2_to_1])
            y <- select_combo[treatment2dose == dose_add, combo_viability]
            curve(
                PharmacoGx::hillCurve(
                    E_ninf = E_ninf_proj,
                    E_inf = E_inf_proj,
                    HS = HS_proj,
                    EC50 = log10(EC50_proj),
                    dose = x
                ),
                from = -10, to = 10, add = TRUE, col = cols[i]
            )
            points(x = log10(dose_to), y = y, col = cols[i])
        }
        if (has_Rsqr[2]) {
            legend(-10, 2,
                legend = paste0("[", treatment2, "] = ", unique_t2_dose,
                                ", R square = ", round(Rsqr_2_to_1, digits = 4)),
                col = cols,
                lty = 1,
                box.lty = 0
            )
        } else {
            legend(-10, 2,
                legend = paste0("[", treatment2, "] = ", unique_t2_dose),
                col = cols,
                lty = 1,
                box.lty = 0
            )
        }
    } else {
        ## unique treatment 1 concentrations
        unique_t1_dose <- unique(select_combo[, treatment1dose])
        cols <- palette(rainbow(length(unique_t1_dose)))
        ## TODO: Find a nicer way to extract R squared value
        if (has_Rsqr[1])
            Rsqr_1_to_2 <- vector(mode = "numeric", length = length(unique_t1_dose))

        ## Initialise an empty background canvas
        plot(
            NULL, xlim = c(-10, 10), ylim = c(0, 2),
            ylab = paste("Response of adding", treatment1, "to", treatment2),
            xlab = paste0("log10([", treatment2,"])"),
            main = title
        )
        for (i in seq_along(unique_t1_dose)) {
            dose_add <- unique_t1_dose[i]
            EC50_proj <- unique(select_combo[treatment1dose == dose_add, EC50_proj_1_to_2])
            HS_proj <- unique(select_combo[treatment1dose == dose_add, HS_proj_1_to_2])
            EC50_add <- unique(select_combo[treatment1dose == dose_add, EC50_1])
            HS_add <- unique(select_combo[treatment1dose == dose_add, HS_1])
            E_inf_add <- unique(select_combo[treatment1dose == dose_add, E_inf_1])
            E_ninf_proj <- PharmacoGx:::.Hill(log10(dose_add), c(HS_add, E_inf_add, log10(EC50_add)))
            E_inf_proj <- unique(select_combo[treatment1dose == dose_add, E_inf_proj_1_to_2])
            dose_to <- select_combo[treatment1dose == dose_add, treatment2dose]
            if (has_Rsqr[1])
                Rsqr_1_to_2[i] <- unique(select_combo[treatment1dose == dose_add, Rsqr_1_to_2])
            y <- select_combo[treatment1dose == dose_add, combo_viability]
            curve(
                PharmacoGx::hillCurve(
                    E_ninf = E_ninf_proj,
                    E_inf = E_inf_proj,
                    HS = HS_proj,
                    EC50 = log10(EC50_proj),
                    dose = x
                ),
                from = -10, to = 10, add = TRUE, col = cols[i]
            )
            points(x = log10(dose_to), y = y, col = cols[i])
        }
        if (has_Rsqr[1]) {
            legend(-10, 2,
                legend = paste0("[", treatment1, "] = ", unique_t1_dose,
                                ", R square = ", round(Rsqr_1_to_2, digits = 4)),
                col = cols,
                lty = 1,
                box.lty = 0
            )
        } else {
            legend(-10, 2,
                legend = paste0("[", treatment1, "] = ", unique_t1_dose),
                col = cols,
                lty = 1,
                box.lty = 0
            )
        }
    }

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
#' @exportMethod computeZIPdelta
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
#' @param residual `character` Method used to minimise residual in fitting curves.
#'     3 methods available: `c("logcosh", "normal", "Cauchy")`.
#'     The default method is `logcosh`.
#'     It minimises the logarithmic hyperbolic cosine loss of the residuals
#'     and provides the fastest estimation among the three methods,
#'     with fitting quality in between `normal` and `Cauchy`;
#'     recommanded when fitting large-scale datasets.
#'     The other two methods minimise residuals by
#'     considering the truncated probability distribution (as in their names) for the residual.
#'     `Cauchy` provides the best fitting quality but also takes the longest to run.
#' @param nthread `integer` Number of cores used to perform computation.
#'     Default 1.
#' @param show_Rsqr `logical` Whether to show the 2-way curve fitting quality in the result.
#'     Default FALSE.
#'
#' @return [TreatmentResponseExperiment] with assay `combo_scores` containing `delta_scores`
#'
#' @examples
#' \dontrun{
#' tre <- computeZIPdelta(tre, residual = "Cauchy", nthread = 2L)
#' }
#'
#' @references
#' Yadav, B., Wennerberg, K., Aittokallio, T., & Tang, J. (2015). Searching for Drug Synergy in Complex Dose–Response Landscapes Using an Interaction Potency Model. Computational and Structural Biotechnology Journal, 13, 504–513. https://doi.org/10.1016/j.csbj.2015.09.001
#'
#' @importFrom CoreGx buildComboProfiles aggregate
#' @importFrom checkmate assertInt assertLogical
#' @import data.table
#' @export
#' @docType methods
setMethod("computeZIPdelta", signature(object = "TreatmentResponseExperiment"),
        function(object, residual = "logcosh", nthread = 1L,
        show_Rsqr = FALSE) {

    if (!is.character(residual)) {
        stop("argument `residual` must be type of logical")
    } else if (length(residual) != 1) {
        stop("argument `residual` must be of length 1")
    }

    assertInt(nthread, lower = 1L)
    assertLogical(show_Rsqr, len = 1)

    combo_keys <- c("treatment1id", "treatment2id",
                    "treatment1dose", "treatment2dose", "sampleid")
    combo_profiles <- tryCatch({
        buildComboProfiles(object, c("HS", "EC50", "E_inf", "ZIP", "combo_viability"))
    }, warning = function(w) {
        message(paste("ZIP reference values have not been pre-computed.",
                      "They will be computed in during delta score calculation."))
        buildComboProfiles(object, c("HS", "EC50", "E_inf", "combo_viability"))
    })
    required_params <- c("HS_1", "HS_2", "E_inf_1", "E_inf_2", "EC50_1", "EC50_2")
    missing_params <- !(required_params %in% colnames(combo_profiles))
    if (any(missing_params))
        stop("Missing required paramters for two-way Hill curve fitting: ",
             paste(required_params[missing_params]))
    has_ZIP <- "ZIP" %in% colnames(combo_profiles)
    if (has_ZIP) {
        combo_ZIP <- combo_profiles[, c(combo_keys, "ZIP"), with = FALSE]
        combo_profiles[, ZIP := NULL]
        setkeyv(combo_ZIP, combo_keys)
    }

    combo_twowayFit <- fitTwowayZIP(combo_profiles = combo_profiles,
                                    residual = residual,
                                    nthread = nthread,
                                    show_Rsqr = show_Rsqr)
    setkeyv(combo_twowayFit, combo_keys)
    if (has_ZIP) {
        combo_twowayFit <- combo_twowayFit[combo_ZIP, on = combo_keys]
        if (show_Rsqr) {
            combo_twowayFit |>
                aggregate(
                    delta_score = .deltaScore(
                        EC50_1_to_2 = EC50_proj_1_to_2,
                        EC50_2_to_1 = EC50_proj_2_to_1,
                        EC50_1 = EC50_1, EC50_2 = EC50_2,
                        HS_1_to_2 = HS_proj_1_to_2,
                        HS_2_to_1 = HS_proj_2_to_1,
                        HS_1 = HS_1, HS_2 = HS_2,
                        E_inf_1 = E_inf_1, E_inf_2 = E_inf_2,
                        E_inf_2_to_1 = E_inf_proj_2_to_1,
                        E_inf_1_to_2 = E_inf_proj_1_to_2,
                        treatment1dose = treatment1dose,
                        treatment2dose = treatment2dose,
                        ZIP = ZIP
                    ),
                    delta_Rsqr_1_to_2 = Rsqr_1_to_2,
                    delta_Rsqr_2_to_1 = Rsqr_2_to_1,
                    by = combo_keys,
                    nthread = nthread
                ) -> delta_scores
        } else {
            combo_twowayFit |>
                aggregate(
                    delta_score = .deltaScore(
                        EC50_1_to_2 = EC50_proj_1_to_2,
                        EC50_2_to_1 = EC50_proj_2_to_1,
                        EC50_1 = EC50_1, EC50_2 = EC50_2,
                        HS_1_to_2 = HS_proj_1_to_2,
                        HS_2_to_1 = HS_proj_2_to_1,
                        HS_1 = HS_1, HS_2 = HS_2,
                        E_inf_1 = E_inf_1, E_inf_2 = E_inf_2,
                        E_inf_2_to_1 = E_inf_proj_2_to_1,
                        E_inf_1_to_2 = E_inf_proj_1_to_2,
                        treatment1dose = treatment1dose,
                        treatment2dose = treatment2dose,
                        ZIP = ZIP
                    ),
                    by = combo_keys,
                    nthread = nthread
                ) -> delta_scores
        }
    } else {
        if (show_Rsqr) {
            combo_twowayFit |>
                aggregate(
                    delta_score = .deltaScore(
                        EC50_1_to_2 = EC50_proj_1_to_2,
                        EC50_2_to_1 = EC50_proj_2_to_1,
                        EC50_1 = EC50_1, EC50_2 = EC50_2,
                        HS_1_to_2 = HS_proj_1_to_2,
                        HS_2_to_1 = HS_proj_2_to_1,
                        HS_1 = HS_1, HS_2 = HS_2,
                        E_inf_1 = E_inf_1, E_inf_2 = E_inf_2,
                        E_inf_2_to_1 = E_inf_proj_2_to_1,
                        E_inf_1_to_2 = E_inf_proj_1_to_2,
                        treatment1dose = treatment1dose,
                        treatment2dose = treatment2dose
                    ),
                    delta_Rsqr_1_to_2 = Rsqr_1_to_2,
                    delta_Rsqr_2_to_1 = Rsqr_2_to_1,
                    by = combo_keys,
                    nthread = nthread
                ) -> delta_scores
        } else {
            combo_twowayFit |>
                aggregate(
                    delta_score = .deltaScore(
                        EC50_1_to_2 = EC50_proj_1_to_2,
                        EC50_2_to_1 = EC50_proj_2_to_1,
                        EC50_1 = EC50_1, EC50_2 = EC50_2,
                        HS_1_to_2 = HS_proj_1_to_2,
                        HS_2_to_1 = HS_proj_2_to_1,
                        HS_1 = HS_1, HS_2 = HS_2,
                        E_inf_1 = E_inf_1, E_inf_2 = E_inf_2,
                        E_inf_2_to_1 = E_inf_proj_2_to_1,
                        E_inf_1_to_2 = E_inf_proj_1_to_2,
                        treatment1dose = treatment1dose,
                        treatment2dose = treatment2dose
                    ),
                    by = combo_keys,
                    nthread = nthread
                ) -> delta_scores
        }
    }
    setkeyv(delta_scores, combo_keys)
    ## Add delta scores to combo_scores in input TRE
    combo_scores <- tryCatch({
        object$combo_scores
    }, error = function(e) {
        NULL
    })
    if (is.null(combo_scores)) {
        ## create a new combo_score assay and save delta scores
        object$combo_scores <- delta_scores
    } else {
        object$combo_scores <- combo_scores[delta_scores, , on = combo_keys]
    }

    return(object)
})

#' @title Vector-based version of [computeZIPdelta]
#'
#' @description
#' Following the calculation of ZIP delta score as in Appendix A3.
#' See reference for details.
#'
#' @param treatment1id `character` a vector of identifiers for treatment 1
#' @param treatment2id `character` a vector of identifiers for treatment 2
#' @param treatment1dose `numeric` a vector of concentrations for treatment 1
#' @param treatment2dose `numeric` a vector of concentrations for treatment 2
#' @param sampleid `character` Cell-line ID of a drug combination screening experiment.
#' @param HS_1 `numeric` Hill coefficient of treatment 1
#' @param HS_2 `numeric` Hill coefficient of treatment 2
#' @param EC50_1 `numeric` relative EC50 of treatment 1.
#' @param EC50_2 `numeric` relative EC50 of treatment 2.
#' @param E_inf_1 `numeric` viability produced by the maximum attainable effect of treatment 1.
#' @param E_inf_2 `numeric` viability produced by the maximum attainable effect of treatment 2.
#' @param combo_viability `numeric` observed viability of the two treatments combined.
#' @param ZIP `numeric` pre-computed ZIP reference values.
#'     If not provided, it will be computed during delta score calculation.
#' @param residual `character` Method used to minimise residual in fitting curves.
#'     3 methods available: `c("logcosh", "normal", "Cauchy")`.
#'     The default method is `logcosh`.
#'     It minimises the logarithmic hyperbolic cosine loss of the residuals
#'     and provides the fastest estimation among the three methods,
#'     with fitting quality in between `normal` and `Cauchy`;
#'     recommanded when fitting large-scale datasets.
#'     The other two methods minimise residuals by
#'     considering the truncated probability distribution (as in their names) for the residual.
#'     `Cauchy` provides the best fitting quality but also takes the longest to run.
#' @param nthread `integer` Number of cores used to perform computation.
#'     Default 1.
#' @param show_Rsqr `logical` Whether to show the 2-way curve fitting quality in the result.
#'     Default FALSE.
#'
#' @return `numeric` delta scores of every dose combinations for any given treatment combinations.
#'
#' @examples
#' \dontrun{
#' ## ZIP is optional. Will be recomputed if not provided.
#' combo_profiles <- CoreGx::buildComboProfiles(tre, c("HS", "EC50", "E_inf", "ZIP", "combo_viability"))
#' combo_profiles[,
#'         .computeZIPdelta(
#'             treatment1id = treatment1id,
#'             treatment2id = treatment2id,
#'             treatment1dose = treatment1dose,
#'             treatment2dose = treatment2dose,
#'             sampleid = sampleid,
#'             HS_1 = HS_1, HS_2 = HS_2,
#'             EC50_1 = EC50_1, EC50_2 = EC50_2,
#'             E_inf_1 = E_inf_1, E_inf_2 = E_inf_2,
#'             combo_viability = combo_viability,
#'             ZIP = ZIP,
#'             nthread = 4,
#'             show_Rsqr = TRUE
#'         )
#'     ] -> delta_scores
#' }
#'
#' @references
#' Yadav, B., Wennerberg, K., Aittokallio, T., & Tang, J. (2015). Searching for Drug Synergy in Complex Dose–Response Landscapes Using an Interaction Potency Model. Computational and Structural Biotechnology Journal, 13, 504–513. https://doi.org/10.1016/j.csbj.2015.09.001
#'
#' @importFrom CoreGx aggregate
#' @importFrom checkmate assertNumeric assertInt assertLogical
#' @import data.table
#' @keywords internal
#' @export
.computeZIPdelta <- function(
    treatment1id, treatment2id, treatment1dose, treatment2dose, sampleid,
    HS_1, HS_2, EC50_1, EC50_2, E_inf_1, E_inf_2, combo_viability, ZIP = NULL,
    residual = "logcosh", nthread = 1L, show_Rsqr = FALSE) {

    assertInt(nthread, lower = 1L)
    assertLogical(show_Rsqr, len = 1)
    len <- length(treatment1dose)
    assertNumeric(treatment1dose, len = len)
    assertNumeric(treatment2dose, len = len)
    assertNumeric(HS_1, len = len)
    assertNumeric(HS_2, len = len)
    assertNumeric(E_inf_1, len = len)
    assertNumeric(E_inf_2, len = len)
    assertNumeric(EC50_1, len = len)
    assertNumeric(EC50_2, len = len)
    assertNumeric(combo_viability, len = len)
    if (!is.null(ZIP))
        assertNumeric(ZIP, len = len)

    combo_keys <- c("treatment1id", "treatment2id",
                    "treatment1dose", "treatment2dose", "sampleid")

    combo_profiles <- data.table(
        treatment1id = treatment1id,
        treatment2id = treatment2id,
        treatment1dose = treatment1dose,
        treatment2dose = treatment2dose,
        sampleid = sampleid,
        combo_viability = combo_viability,
        HS_1 = HS_1,
        HS_2 = HS_2,
        EC50_1 = EC50_1,
        EC50_2 = EC50_2,
        E_inf_1 = E_inf_1,
        E_inf_2 = E_inf_2
    )

    if (!is.null(ZIP)) {
        combo_ZIP <- data.table(
            treatment1id = treatment1id,
            treatment2id = treatment2id,
            treatment1dose = treatment1dose,
            treatment2dose = treatment2dose,
            sampleid = sampleid,
            ZIP = ZIP
        )
        setkeyv(combo_ZIP, combo_keys)
    }

    combo_twowayFit <- fitTwowayZIP(combo_profiles = combo_profiles,
                                    residual = residual,
                                    nthread = nthread,
                                    show_Rsqr = show_Rsqr)
    setkeyv(combo_twowayFit, combo_keys)
    if (is.null(ZIP)) {
        if (show_Rsqr) {
            combo_twowayFit |>
                aggregate(
                    delta_score = .deltaScore(
                        EC50_1_to_2 = EC50_proj_1_to_2,
                        EC50_2_to_1 = EC50_proj_2_to_1,
                        EC50_1 = EC50_1, EC50_2 = EC50_2,
                        HS_1_to_2 = HS_proj_1_to_2,
                        HS_2_to_1 = HS_proj_2_to_1,
                        HS_1 = HS_1, HS_2 = HS_2,
                        E_inf_1 = E_inf_1, E_inf_2 = E_inf_2,
                        E_inf_2_to_1 = E_inf_proj_2_to_1,
                        E_inf_1_to_2 = E_inf_proj_1_to_2,
                        treatment1dose = treatment1dose,
                        treatment2dose = treatment2dose
                    ),
                    delta_Rsqr_1_to_2 = Rsqr_1_to_2,
                    delta_Rsqr_2_to_1 = Rsqr_2_to_1,
                    by = combo_keys,
                    nthread = nthread
                ) -> delta_scores
        } else {
            combo_twowayFit |>
                aggregate(
                    delta_score = .deltaScore(
                        EC50_1_to_2 = EC50_proj_1_to_2,
                        EC50_2_to_1 = EC50_proj_2_to_1,
                        EC50_1 = EC50_1, EC50_2 = EC50_2,
                        HS_1_to_2 = HS_proj_1_to_2,
                        HS_2_to_1 = HS_proj_2_to_1,
                        HS_1 = HS_1, HS_2 = HS_2,
                        E_inf_1 = E_inf_1, E_inf_2 = E_inf_2,
                        E_inf_2_to_1 = E_inf_proj_2_to_1,
                        E_inf_1_to_2 = E_inf_proj_1_to_2,
                        treatment1dose = treatment1dose,
                        treatment2dose = treatment2dose
                    ),
                    by = combo_keys,
                    nthread = nthread
                ) -> delta_scores
        }
    } else {
        combo_twowayFit <- combo_twowayFit[combo_ZIP, on = combo_keys]
        if (show_Rsqr) {
            combo_twowayFit |>
                aggregate(
                    delta_score = .deltaScore(
                        EC50_1_to_2 = EC50_proj_1_to_2,
                        EC50_2_to_1 = EC50_proj_2_to_1,
                        EC50_1 = EC50_1, EC50_2 = EC50_2,
                        HS_1_to_2 = HS_proj_1_to_2,
                        HS_2_to_1 = HS_proj_2_to_1,
                        HS_1 = HS_1, HS_2 = HS_2,
                        E_inf_1 = E_inf_1, E_inf_2 = E_inf_2,
                        E_inf_2_to_1 = E_inf_proj_2_to_1,
                        E_inf_1_to_2 = E_inf_proj_1_to_2,
                        treatment1dose = treatment1dose,
                        treatment2dose = treatment2dose,
                        ZIP = ZIP
                    ),
                    delta_Rsqr_1_to_2 = Rsqr_1_to_2,
                    delta_Rsqr_2_to_1 = Rsqr_2_to_1,
                    by = combo_keys,
                    nthread = nthread
                ) -> delta_scores
        } else {
            combo_twowayFit |>
                aggregate(
                    delta_score = .deltaScore(
                        EC50_1_to_2 = EC50_proj_1_to_2,
                        EC50_2_to_1 = EC50_proj_2_to_1,
                        EC50_1 = EC50_1, EC50_2 = EC50_2,
                        HS_1_to_2 = HS_proj_1_to_2,
                        HS_2_to_1 = HS_proj_2_to_1,
                        HS_1 = HS_1, HS_2 = HS_2,
                        E_inf_1 = E_inf_1, E_inf_2 = E_inf_2,
                        E_inf_2_to_1 = E_inf_proj_2_to_1,
                        E_inf_1_to_2 = E_inf_proj_1_to_2,
                        treatment1dose = treatment1dose,
                        treatment2dose = treatment2dose,
                        ZIP = ZIP
                    ),
                    by = combo_keys,
                    nthread = nthread
                ) -> delta_scores
        }
    }
    if (show_Rsqr) {
        return(as.list(delta_scores[,
            c(combo_keys, "delta_score", "delta_Rsqr_1_to_2", "delta_Rsqr_2_to_1"),
            with = FALSE
        ]))
    } else {
        return(as.list(delta_scores[, c(combo_keys, "delta_score"), with = FALSE]))
    }
}


#' @title Calculate ZIP delta score for a drug combination
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
#' @param E_inf_2_to_1 `numeric` projected maximum attainable effect of
#'     adding treatment 2 to treatment 1.
#' @param E_inf_1_to_2 `numeric` projected maximum attainable effect of
#'     adding treatment 1 to treatment 2.
#' @param E_inf_1 `numeric` viability produced by the maximum attainable effect of treatment 1.
#' @param E_inf_2 `numeric` viability produced by the maximum attainable effect of treatment 2.
#' @param treatment1dose `numeric` a vector of concentrations for treatment 1
#' @param treatment2dose `numeric` a vector of concentrations for treatment 2
#' @param ZIP `numeric` pre-computed ZIP reference values.
#'     If not provided, it will be computed during delta score calculation.
#'
#' @return `numeric` a ZIP delta score to quantify synergy for the drug combination.
#'
#' @noRd
#' @references
#' Yadav, B., Wennerberg, K., Aittokallio, T., & Tang, J. (2015). Searching for Drug Synergy in Complex Dose–Response Landscapes Using an Interaction Potency Model. Computational and Structural Biotechnology Journal, 13, 504–513. https://doi.org/10.1016/j.csbj.2015.09.001
#' @export
.deltaScore <- function(EC50_1_to_2, EC50_2_to_1, EC50_1, EC50_2,
                        HS_1_to_2, HS_2_to_1, HS_1, HS_2,
                        E_inf_2_to_1, E_inf_1_to_2, E_inf_1, E_inf_2,
                        treatment1dose, treatment2dose, ZIP = NULL) {

    viability_1 <- .Hill(log10(treatment1dose), c(HS_1, E_inf_1, log10(EC50_1)))
    viability_2 <- .Hill(log10(treatment2dose), c(HS_2, E_inf_2, log10(EC50_2)))
    viability_2_to_1 <- hillCurve(
        dose = treatment1dose,
        HS = HS_2_to_1,
        EC50 = EC50_2_to_1,
        E_ninf = viability_2,
        E_inf = E_inf_2_to_1
    )
    viability_1_to_2 <- hillCurve(
        dose = treatment2dose,
        HS = HS_1_to_2,
        EC50 = EC50_1_to_2,
        E_ninf = viability_1,
        E_inf = E_inf_1_to_2
    )
    ## avoid re-calculating ZIP references
    if (is.null(ZIP)) {
        viability_ZIP <- computeZIP(treatment1dose = treatment1dose,
                                    treatment2dose = treatment2dose,
                                    HS_1 = HS_1, HS_2 = HS_2,
                                    EC50_1 = EC50_1, EC50_2 = EC50_2,
                                    E_inf_1 = E_inf_1, E_inf_2 = E_inf_2)
    } else {
        viability_ZIP <- ZIP
    }
    delta <- viability_ZIP - (1/2) * (viability_2_to_1 + viability_1_to_2)

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
#' @examples
#' (bliss <- computeBliss(0.75, 0.65))
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
#' @examples
#' (hsa <- computeHSA(0.75, 0.65))
#'
#' @export
computeHSA <- function(viability_1, viability_2) {
    HSA_ref <- min(viability_1, viability_2)
    return(HSA_ref)
}
