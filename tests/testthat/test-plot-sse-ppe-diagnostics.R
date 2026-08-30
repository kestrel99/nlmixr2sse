skip_on_cran()

# ---- Cramer-von Mises discrepancy -------------------------------------

test_that("Cramer-von Mises separates a correct model from a misspecified one (power)", {
  good <- ppe_dofv(n = 400L, df = 1, ncp = 9, seed = 31L)
  bad <- c(ppe_dofv(200L, df = 1, ncp = 2, seed = 32L),
           ppe_dofv(200L, df = 1, ncp = 30, seed = 33L))

  cvm_good <- .ppeCramerVonMises(good[good > 0], df = 1, ncp = 9)
  cvm_bad <- .ppeCramerVonMises(bad[bad > 0], df = 1,
                                ncp = .ppeChiSquareMle(bad, df = 1)$estimate)

  expect_lt(cvm_good, cvm_bad)
})

test_that("Cramer-von Mises separates a correct model from a misspecified one (type1)", {
  # target = "df": ncp fixed at 0, df is the value being diagnosed -- the
  # Type-I-comparison counterpart to the power-mode test above. "Correct"
  # here is a single central chi-square(df = 3); "misspecified" is a mixture
  # of two central chi-squares with very different df, which a single fitted
  # df cannot describe well.
  good <- ppe_dofv(n = 400L, df = 3, ncp = 0, seed = 51L)
  bad <- c(ppe_dofv(200L, df = 1, ncp = 0, seed = 53L),
           ppe_dofv(200L, df = 10, ncp = 0, seed = 54L))

  fit_good <- .ppeChiSquareMle(good, ncp = 0)
  fit_bad <- .ppeChiSquareMle(bad, ncp = 0)

  cvm_good <- .ppeCramerVonMises(good[good > 0], df = fit_good$estimate, ncp = 0, target = "df")
  cvm_bad <- .ppeCramerVonMises(bad[bad > 0], df = fit_bad$estimate, ncp = 0, target = "df")

  expect_lt(cvm_good, cvm_bad)
})

# ---- Diagnostic p-value -------------------------------------------------

test_that("the diagnostic p-value is small for a misspecified mixture (power)", {
  bad <- c(ppe_dofv(200L, df = 1, ncp = 2, seed = 34L),
           ppe_dofv(200L, df = 1, ncp = 30, seed = 35L))
  fit <- .ppeChiSquareMle(bad, df = 1)

  p <- .ppeDiagnosticPValue(bad[bad > 0], df = 1, estimate = fit$estimate,
                            bootstrapSamples = 200L, seed = 36L)

  expect_lt(p, 0.05)
})

test_that("the diagnostic p-value is not small when the model is correct (power)", {
  good <- ppe_dofv(n = 400L, df = 1, ncp = 9, seed = 37L)
  fit <- .ppeChiSquareMle(good, df = 1)

  p <- .ppeDiagnosticPValue(good[good > 0], df = 1, estimate = fit$estimate,
                            bootstrapSamples = 200L, seed = 38L)

  expect_gt(p, 0.05)
})

test_that("the diagnostic p-value is small for a misspecified mixture (type1)", {
  bad <- c(ppe_dofv(200L, df = 1, ncp = 0, seed = 53L),
           ppe_dofv(200L, df = 10, ncp = 0, seed = 54L))
  fit <- .ppeChiSquareMle(bad, ncp = 0)

  p <- .ppeDiagnosticPValue(bad[bad > 0], df = NULL, estimate = fit$estimate,
                            target = "df", bootstrapSamples = 200L, seed = 55L)

  expect_lt(p, 0.05)
})

test_that("the diagnostic p-value is not small when the model is correct (type1)", {
  good <- ppe_dofv(n = 400L, df = 3, ncp = 0, seed = 51L)
  fit <- .ppeChiSquareMle(good, ncp = 0)

  p <- .ppeDiagnosticPValue(good[good > 0], df = NULL, estimate = fit$estimate,
                            target = "df", bootstrapSamples = 200L, seed = 52L)

  expect_gt(p, 0.05)
})

# ---- plotSSEPpeDiagnostics() ---------------------------------------------

test_that("plotSSEPpeDiagnostics returns a ggplot carrying auditable data", {
  skip_if_not_installed("patchwork")
  sse <- fake_ppe_sse_object(df = 1, ncp = 10, n = 80L, seed = 41L)
  p <- plotSSEPpeDiagnostics(
    sse, comparisons = sseComparison("simulation", "alt1", df = 1),
    bootstrapSamples = 50L, bootSeed = 2L
  )

  expect_s3_class(p, "ggplot")
  diag <- attr(p, "ppeDiagnostics")
  expect_true(all(c("comparison", "cvm", "p_value", "n", "n_nonpositive") %in%
                    names(diag)))
  expect_equal(nrow(diag), 1L)
})

test_that("the diagnostic subtitle names the excluded count", {
  skip_if_not_installed("patchwork")
  sse <- fake_ppe_sse_object(df = 1, ncp = 10, n = 80L, seed = 42L,
                             nNonPositive = 6L)
  p <- suppressWarnings(plotSSEPpeDiagnostics(
    sse, comparisons = sseComparison("simulation", "alt1", df = 1),
    bootstrapSamples = 20L, bootSeed = 2L
  ))

  expect_match(p$labels$subtitle, "6")
})

test_that("the caption states the envelope is pointwise, not simultaneous", {
  skip_if_not_installed("patchwork")
  sse <- fake_ppe_sse_object(df = 1, ncp = 10, n = 60L, seed = 43L)
  p <- plotSSEPpeDiagnostics(
    sse, comparisons = sseComparison("simulation", "alt1", df = 1),
    bootstrapSamples = 20L, bootSeed = 2L
  )

  expect_match(p$labels$caption, "pointwise")
  expect_match(p$labels$caption, "not simultaneous")
})

test_that("multiple comparisons produce one row per comparison, distinguishing power from type1", {
  skip_if_not_installed("patchwork")
  power_sse <- fake_ppe_sse_object(df = 1, ncp = 10, n = 60L, seed = 44L)
  type1_sse <- fake_ppe_type1_sse_object(df = 1, n = 60L, seed = 44L)

  power_cmp <- sseComparison("simulation", "alt1", df = 1, label = "power check")
  type1_cmp <- sseComparison("alt1", "simulation", df = 1, label = "type1 check")

  p_power <- plotSSEPpeDiagnostics(power_sse, comparisons = power_cmp,
                                   bootstrapSamples = 20L, bootSeed = 5L)
  p_type1 <- plotSSEPpeDiagnostics(type1_sse, comparisons = type1_cmp,
                                   bootstrapSamples = 20L, bootSeed = 5L)

  diag <- rbind(attr(p_power, "ppeDiagnostics"), attr(p_type1, "ppeDiagnostics"))
  expect_equal(nrow(diag), 2L)
  expect_equal(diag$comparison, c("power check", "type1 check"))

  # And a single call with both comparisons resolved together also produces
  # one row per comparison, faceted rather than collapsed.
  p_both <- plotSSEPpeDiagnostics(
    power_sse, comparisons = list(power_cmp), bootstrapSamples = 20L, bootSeed = 5L
  )
  expect_equal(nrow(attr(p_both, "ppeDiagnostics")), 1L)
})

test_that("a criticalValue-only comparison is refused", {
  skip_if_not_installed("patchwork")
  sse <- fake_ppe_sse_object(df = 1, ncp = 10, n = 40L, seed = 45L)
  cmp <- sseComparison("simulation", "alt1", criticalValue = 5)

  expect_error(plotSSEPpeDiagnostics(sse, comparisons = cmp), "noncentral chi-square")
})
