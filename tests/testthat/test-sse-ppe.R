skip_on_cran()

test_that("the MLE recovers a known noncentrality parameter", {
  for (spec in list(list(df = 1, ncp = 8), list(df = 3, ncp = 15), list(df = 5, ncp = 25))) {
    d <- ppe_dofv(n = 4000L, df = spec$df, ncp = spec$ncp, seed = 7L)
    fit <- .ppeChiSquareMle(d, df = spec$df)
    expect_equal(fit$estimate, spec$ncp, tolerance = 0.05)
  }
})

test_that("the MLE recovers a known df when ncp is held at zero", {
  d <- ppe_dofv(n = 4000L, df = 4, ncp = 0, seed = 11L)
  fit <- .ppeChiSquareMle(d, ncp = 0)
  expect_equal(fit$estimate, 4, tolerance = 0.05)
})

test_that("the MLE agrees with a direct likelihood grid", {
  d <- ppe_dofv(n = 300L, df = 2, ncp = 6, seed = 13L)
  grid <- seq(0.1, 20, by = 0.01)
  ll <- vapply(grid, function(ncp) sum(stats::dchisq(d[d > 0], 2, ncp, log = TRUE)), numeric(1))

  expect_equal(.ppeChiSquareMle(d, df = 2)$estimate, grid[which.max(ll)], tolerance = 0.02)
})

test_that("non-positive values are counted, never hidden", {
  d <- c(-2, -1, 0, 3, 4, 5, 6)
  fit <- .ppeChiSquareMle(d, df = 1)

  expect_equal(fit$nRetained, 4L)
  expect_equal(fit$nNonPositive, 3L)
  expect_equal(fit$n, 7L)
})

test_that("a boundary solution is flagged rather than reported as a clean fit", {
  d <- c(0.1, 0.2, 0.15, 0.3)
  fit <- .ppeChiSquareMle(d, df = 4)

  expect_true(fit$boundary)
  expect_lt(fit$estimate, 1e-6)
})

test_that("too few positive values is an error, not a silent fit", {
  expect_error(.ppeChiSquareMle(c(-1, 2), df = 1), "at least 2")
})

test_that(".ppeChiSquareMle rejects both df and ncp supplied together", {
  expect_error(.ppeChiSquareMle(c(1, 2, 3), df = 1, ncp = 2), "exactly one")
})

test_that(".ppeChiSquareMle rejects neither df nor ncp supplied", {
  expect_error(.ppeChiSquareMle(c(1, 2, 3)), "exactly one")
})

test_that("the nonpositive policy controls warning and abort, never the counts", {
  expect_warning(.ppeApplyNonpositivePolicy(2L, 5L, "warn"), "2 of 5")
  expect_error(.ppeApplyNonpositivePolicy(2L, 5L, "error"), "2 of 5")
  expect_silent(.ppeApplyNonpositivePolicy(2L, 5L, "drop"))
})

# --- .ppePowerPlotData(method = "distribution_mle") ---------------------

test_that("the fitted noncentrality does not depend on the threshold, unlike exceedance", {
  # This is the defining difference from method = "exceedance": one
  # noncentrality fitted from the whole distribution, reused for every
  # threshold, rather than a fresh implied noncentrality per threshold.
  #
  # df = 2 is passed explicitly via sseComparison(): fake_ppe_sse_object()'s
  # fixture model spec (fake_sse_fit_specs()) always implies df = 1 from its
  # own parameter counts, independent of the `df` used to generate the draws
  # -- so recovering the true generating ncp needs the true df supplied too.
  sse <- fake_ppe_sse_object(df = 2, ncp = 10, n = 500L, seed = 21L)
  cmp <- sseComparison("simulation", "alt1", df = 2)
  data <- .ppePowerPlotData(
    sse,
    thresholds = c(4, 12),
    studySizes = 12L,
    comparisons = cmp,
    nonpositivePolicy = "drop"
  )

  ncp_by_threshold <- split(data$mle_ncp, data$threshold)
  expect_length(unique(vapply(ncp_by_threshold, unique, numeric(1))), 1L)

  # And, unlike exceedance, that single fitted value should recover the
  # generating ncp (loosely -- this is a sanity check, not a tolerance test;
  # the tight tolerance checks already live in the .ppeChiSquareMle() tests).
  expect_equal(unique(data$mle_ncp), 10, tolerance = 0.1)
})

test_that("distribution_mle refuses a comparison with an explicit criticalValue", {
  sse <- fake_ppe_sse_object(df = 1, ncp = 8, n = 20L, seed = 5L)
  cmp <- sseComparison("simulation", "alt1", criticalValue = 5)

  expect_error(
    .ppePowerPlotData(sse, comparisons = cmp, nonpositivePolicy = "drop"),
    "not eligible"
  )
})

test_that("distribution_mle warns when df is inferred rather than declared", {
  # distribution_mle treats df as known -- it is the exponent in the
  # noncentral chi-square density the whole fit rests on -- so, unlike the
  # legacy exceedance branch, an inferred df must not be silent.
  sse <- fake_ppe_sse_object(df = 1, ncp = 8, n = 20L, seed = 5L)

  expect_warning(
    .ppePowerPlotData(sse, nonpositivePolicy = "drop"),
    "Degrees of freedom inferred"
  )
})

test_that("df_source reports explicit vs. inferred provenance", {
  sse <- fake_ppe_sse_object(df = 1, ncp = 8, n = 20L, seed = 5L)

  explicit_cmp <- sseComparison("simulation", "alt1", df = 1)
  explicit_data <- .ppePowerPlotData(
    sse,
    comparisons = explicit_cmp,
    nonpositivePolicy = "drop"
  )
  expect_true(all(explicit_data$df_source == "explicit"))

  inferred_data <- suppressWarnings(
    .ppePowerPlotData(sse, nonpositivePolicy = "drop")
  )
  expect_true(all(inferred_data$df_source == "parameter_count"))
})
