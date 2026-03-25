library(PharmacoGx)

context("Testing drugDoseResponseCurve.")

capture_curve_calls <- function() {
  calls <- new.env(parent = emptyenv())
  calls$points <- list()
  calls$lines <- list()
  calls$fits <- list()

  list(
    calls = calls,
    points = function(x, y, ...) {
      calls$points[[length(calls$points) + 1]] <- list(
        x = x,
        y = y,
        args = list(...)
      )
      invisible(NULL)
    },
    lines = function(x, y, ...) {
      calls$lines[[length(calls$lines) + 1]] <- list(
        x = x,
        y = y,
        args = list(...)
      )
      invisible(NULL)
    },
    fit = function(conc, viability, ...) {
      calls$fits[[length(calls$fits) + 1]] <- list(
        conc = conc,
        viability = viability,
        args = list(...)
      )
      list(HS = 1, E0 = 100, E_inf = 0, EC50 = 1)
    }
  )
}

open_test_device <- function() {
  grDevices::pdf(file = tempfile(fileext = ".pdf"))
}

make_identical_grid_fixture <- function() {
  data("GDSCsmall", package = "PharmacoGx")
  pSet <- subsetTo(GDSCsmall, cells = "22RV1", drugs = "AZD6482")
  stopifnot(nrow(sensitivityInfo(pSet)) == 2)

  raw <- sensitivityRaw(pSet)
  raw[2, , "Dose"] <- raw[1, , "Dose"]
  raw[1, , "Viability"] <- as.character(c(100, 90, 80, 70, 60, 50, 40, 30, 20))
  raw[2, , "Viability"] <- as.character(c(80, 70, 60, 50, 40, 30, 20, 10, 0))
  sensitivityRaw(pSet) <- raw

  pSet
}

test_that("summarize.replicates collapses duplicate doses before plotting and fitting", {
  pSet <- make_identical_grid_fixture()
  expected.dose <- as.numeric(sensitivityRaw(pSet)[1, , "Dose"])
  expected.viability <- c(90, 80, 70, 60, 50, 40, 30, 20, 10)
  capture <- capture_curve_calls()

  open_test_device()
  on.exit(grDevices::dev.off(), add = TRUE)

  testthat::local_mocked_bindings(
    points = capture$points,
    lines = capture$lines,
    logLogisticRegression = capture$fit,
    .package = "PharmacoGx"
  )

  expect_silent(
    expect_invisible(
      drugDoseResponseCurve(
        drug = "AZD6482",
        cellline = "22RV1",
        pSets = pSet,
        plot.type = "Both",
        summarize.replicates = TRUE
      )
    )
  )

  expect_length(capture$calls$points, 1)
  expect_equal(as.numeric(capture$calls$points[[1]]$x), expected.dose)
  expect_equal(as.numeric(capture$calls$points[[1]]$y), expected.viability)

  expect_length(capture$calls$fits, 1)
  expect_equal(as.numeric(capture$calls$fits[[1]]$conc), expected.dose)
  expect_equal(
    as.numeric(capture$calls$fits[[1]]$viability),
    expected.viability
  )

  expect_gte(length(capture$calls$lines), 2)
  expect_equal(as.numeric(capture$calls$lines[[1]]$x), expected.dose)
  expect_equal(as.numeric(capture$calls$lines[[1]]$y), expected.viability)
})

test_that("summarize.replicates warns when replicate dose grids differ", {
  data("GDSCsmall", package = "PharmacoGx")
  capture <- capture_curve_calls()

  open_test_device()
  on.exit(grDevices::dev.off(), add = TRUE)

  testthat::local_mocked_bindings(
    points = capture$points,
    lines = capture$lines,
    .package = "PharmacoGx"
  )

  expect_warning(
    expect_invisible(
      drugDoseResponseCurve(
        drug = "AZD6482",
        cellline = "22RV1",
        pSets = GDSCsmall,
        plot.type = "Actual",
        summarize.replicates = TRUE
      )
    ),
    regexp = paste(
      "Replicate dose grids differ for AZD6482:22RV1.*",
      "only exact duplicate doses were summarized."
    )
  )

  expect_length(capture$calls$points, 1)
  expect_length(capture$calls$points[[1]]$x, 18)
})

test_that("summarize.replicates = FALSE leaves replicate observations unsummarized", {
  pSet <- make_identical_grid_fixture()
  expected.dose <- as.numeric(sensitivityRaw(pSet)[1, , "Dose"])
  first.viability <- c(100, 90, 80, 70, 60, 50, 40, 30, 20)
  second.viability <- c(80, 70, 60, 50, 40, 30, 20, 10, 0)
  capture <- capture_curve_calls()

  open_test_device()
  on.exit(grDevices::dev.off(), add = TRUE)

  testthat::local_mocked_bindings(
    points = capture$points,
    lines = capture$lines,
    .package = "PharmacoGx"
  )

  expect_silent(
    expect_invisible(
      drugDoseResponseCurve(
        drug = "AZD6482",
        cellline = "22RV1",
        pSets = pSet,
        plot.type = "Actual",
        summarize.replicates = FALSE
      )
    )
  )

  expect_length(capture$calls$points, 2)
  expect_equal(as.numeric(capture$calls$points[[1]]$x), expected.dose)
  expect_equal(as.numeric(capture$calls$points[[1]]$y), first.viability)
  expect_equal(as.numeric(capture$calls$points[[2]]$x), expected.dose)
  expect_equal(as.numeric(capture$calls$points[[2]]$y), second.viability)

  expect_length(capture$calls$lines, 2)
})

test_that("manual concentrations and viabilities plot with default legend labels", {
  open_test_device()
  on.exit(grDevices::dev.off(), add = TRUE)

  expect_silent(
    expect_invisible(
      drugDoseResponseCurve(
        concentrations = list("Experiment 1" = c(0.008, 0.04, 0.2, 1)),
        viabilities = list(c(100, 50, 30, 1)),
        plot.type = "Both"
      )
    )
  )
})

test_that("manual unsorted concentrations warn but still plot", {
  open_test_device()
  on.exit(grDevices::dev.off(), add = TRUE)

  expect_warning(
    expect_invisible(
      drugDoseResponseCurve(
        concentrations = list("Experiment 1" = c(1, 0.2, 0.04, 0.008)),
        viabilities = list(c(1, 30, 50, 100)),
        plot.type = "Actual"
      )
    ),
    regexp = "Concentration Values were unsorted"
  )
})
