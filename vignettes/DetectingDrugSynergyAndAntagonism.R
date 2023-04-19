## ----eval=TRUE, include=FALSE-------------------------------------------------
# convenience variables
cgx <- BiocStyle::Biocpkg("CoreGx")
pgx <- BiocStyle::Biocpkg("PharmacoGx")
dt <- BiocStyle::CRANpkg("data.table")

# knitr options
knitr::opts_chunk$set(warning=FALSE)

## ----load_dependencies_eval, eval=TRUE, echo=FALSE----------------------------
suppressPackageStartupMessages({
    library(PharmacoGx)
    library(CoreGx)
    library(data.table)
    library(ggplot2)
})

## ----load_dependencies_echo, eval=FALSE, echo=TRUE----------------------------
#  library(PharmacoGx)
#  library(CoreGx)
#  library(data.table)
#  library(ggplot2)

## -----------------------------------------------------------------------------
input_file <- system.file("extdata/mathews_griner.csv.tar.gz",
    package="PharmacoGx")
mathews_griner <- fread(input_file)

## ----experiment_design_hypothesis---------------------------------------------
groups <- list(
    rowDataMap=c(
        treatment1id="RowName", treatment2id="ColName",
        treatment1dose="RowConcs", treatment2dose="ColConcs"
    ),
    colDataMap=c("sampleid")
)
groups[["assayMap"]] <- c(groups$rowDataMap, groups$colDataMap)
(groups)

## ----handling_technical_replicates--------------------------------------------
# The := operator modifies a data.table by reference (i.e., without making a copy)
mathews_griner[, tech_rep := seq_len(.N), by=c(groups[["assayMap"]])]
if (max(mathews_griner[["tech_rep"]]) > 1) {
    groups[["colDataMap"]] <- c(groups[["colDataMap"]], "tech_rep")
    groups[["assayMap"]] <- c(groups[["assayMap"]], "tech_rep")
} else {
    # delete the additional column if not needed
    message("No technical replicates in this dataset!")
    mathews_griner[["tech_reps"]] <- NULL
}

## ----build_tredatamapper------------------------------------------------------
(treMapper <- TREDataMapper(rawdata=mathews_griner))

## ----evaluate_tre_mapping_guess-----------------------------------------------
(guess <- guessMapping(treMapper, groups, subset=TRUE))

## ----update_tredatamapper_with_guess------------------------------------------
metadataMap(treMapper) <- list(experiment_metadata=guess$metadata$mapped_columns)
rowDataMap(treMapper) <- guess$rowDataMap
colDataMap(treMapper) <- guess$colDataMap
assayMap(treMapper) <- list(raw=guess$assayMap)
treMapper

## ----metaconstruct_the_tre----------------------------------------------------
(tre <- metaConstruct(treMapper))

## ----normalize_to_dose_0_0_control--------------------------------------------
raw <- tre[["raw"]]
raw[,
    viability := viability / .SD[treatment1dose == 0 & treatment2dose == 0, viability],
    by=c("treatment1id", "treatment2id", "sampleid", "tech_rep")
]
raw[, viability := pmax(0, viability)]  # truncate min viability at 0
tre[["raw"]] <- raw

## ----sanity_check_viability---------------------------------------------------
tre[["raw"]][, range(viability)]

## ----find_bad_viability_treatment, warning=FALSE------------------------------
(bad_treatments <- tre[["raw"]][viability > 2, unique(treatment1id)])

## ----remove_bad_viability_treatment, warning=FALSE----------------------------
(tre <- subset(tre, !(treatment1id %in% bad_treatments)))

## ----sanity_check_viability2--------------------------------------------------
tre[["raw"]][, range(viability)]

## ----creating_monotherapy_assay-----------------------------------------------
tre_qc <- tre |>
    endoaggregate(
        subset=treatment2dose == 0,  # filter to only monotherapy rows
        assay="raw",
        target="mono_viability",  # create a new assay named mono_viability
        mean_viability=pmin(1, mean(viability)),
        by=c("treatment1id", "treatment1dose", "sampleid")
    )

## ----monotherapy_curve_fits, messages=FALSE-----------------------------------
tre_fit <- tre_qc |>
    endoaggregate(
        {  # the entire code block is evaluated for each group in our group by
            # 1. fit a log logistic curve over the dose range
            fit <- PharmacoGx::logLogisticRegression(treatment1dose, mean_viability,
                viability_as_pct=FALSE)
            # 2. compute curve summary metrics
            ic50 <- computeIC50(treatment1dose, Hill_fit=fit)
            aac <- computeAUC(treatment1dose, Hill_fit=fit)
            # 3. assemble the results into a list, each item will become a
            #   column in the target assay.
            list(
                HS=fit[["HS"]],
                E_inf = fit[["E_inf"]],
                EC50 = fit[["EC50"]],
                Rsq=as.numeric(unlist(attributes(fit))),
                aac_recomputed=aac,
                ic50_recomputed=ic50
            )
        },
        assay="mono_viability",
        target="mono_profiles",
        enlist=FALSE,  # this option enables the use of a code block for aggregation
        by=c("treatment1id", "sampleid"),
        nthread=2  # parallelize over multiple cores to speed up the computation
    )

## ----create_combo_viability, message=FALSE------------------------------------
tre_combo <- tre_fit |>
    endoaggregate(
        assay="raw",
        target="combo_viability",
        mean(viability),
        by=c("treatment1id", "treatment2id", "treatment1dose", "treatment2dose",
            "sampleid")
    )

## ----add_monotherapy_fits_to_combo_viability----------------------------------
tre_combo <- tre_combo |>
    mergeAssays(
        x="combo_viability",
        y="mono_profiles",
        by=c("treatment1id", "sampleid")
    ) |>
    mergeAssays(
        x="combo_viability",
        y="mono_profiles",
        by.x=c("treatment2id", "sampleid"),
        by.y=c("treatment1id", "sampleid"),
        suffixes=c("_1", "_2")  # add sufixes to duplicate column names
    )

## -----------------------------------------------------------------------------
tre_combo <- tre_combo |>
    endoaggregate(
        viability_1=.SD[treatment2dose == 0, mean_viability],
        assay="combo_viability",
        by=c("treatment1id", "treatment1dose", "sampleid")
    ) |>
    endoaggregate(
        viability_2=.SD[treatment1dose == 0, mean_viability],
        assay="combo_viability",
        by=c("treatment1id", "treatment2dose", "sampleid")
    )

## ----compute_synergy_null_hypotheses, message=FALSE---------------------------
tre_synergy <- tre_combo |>
    endoaggregate(
        assay="combo_viability",
        HSA_ref=computeHSA(viability_1, viability_2),
        Bliss_ref=computeBliss(viability_1, viability_2),
        Loewe_ref=computeLoewe(
            treatment1dose, HS_1=HS_1, EC50_1=EC50_1, E_inf_1=E_inf_1,
            treatment2dose, HS_2=HS_2, EC50_2=EC50_2, E_inf_2=E_inf_2
        ),
        ZIP_ref=computeZIP(
            treatment1dose, HS_1=HS_1, EC50_1=EC50_1, E_inf_1=E_inf_1,
            treatment2dose, HS_2=HS_2, EC50_2=EC50_2, E_inf_2=E_inf_2
        ),
        by=assayKeys(tre_combo, "combo_viability"),
        nthread=2
    )

## ----synergy_score_vs_reference-----------------------------------------------
tre_synergy <- tre_synergy |>
    endoaggregate(
        assay="combo_viability",
        HSA_score=HSA_ref - mean_viability,
        Bliss_score=Bliss_ref - mean_viability,
        Loewe_score=Loewe_ref - mean_viability,
        ZIP_score=ZIP_ref - mean_viability,
        by=assayKeys(tre_synergy, "combo_viability")
    )

