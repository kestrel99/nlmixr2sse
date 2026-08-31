skip_on_cran()

# ---- RNG isolation --------------------------------------------------------

test_that(".ppeEnvelopeBounds leaves the caller's RNG stream untouched", {
  # Mirrors the existing .ppeParametricBootstrap() isolation test: the
  # envelope draws its own bootstrap sample internally and must restore
  # .Random.seed (including its absence) exactly as .withPpeSeed() does
  # elsewhere in the PPE code, not just "usually leave it alone".
  set.seed(99L)
  before <- .Random.seed

  .ppeEnvelopeBounds(
    n = 50L, df = 1, ncp = 8, grid = seq(0, 20, length.out = 50L),
    envelopeSamples = 20L, conf.level = 0.95, seed = 5L
  )

  expect_equal(.Random.seed, before)
})

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

test_that("the diagnostic p-value cannot be exactly zero", {
  # mean(null_stats >= observed) can be exactly 0 with a finite bootstrap
  # sample, which is not a valid Monte Carlo p-value (it implies infinite
  # certainty against the fitted model from a finite simulation). The
  # standard correction is (1 + #{T_b >= T_obs}) / (B + 1), which is always
  # strictly positive.
  good <- ppe_dofv(n = 200L, df = 4, ncp = 15, seed = 81L)
  fit <- .ppeChiSquareMle(good, df = 4)

  p <- .ppeDiagnosticPValue(
    fit$retained, df = 4, estimate = fit$estimate, target = "ncp",
    bootstrapSamples = 20L, seed = 5L
  )

  expect_gt(p, 0)
  expect_true(p >= 1 / 21)
})

test_that("the diagnostic p-value floor is 1/(bootstrapSamples + 1), never 0", {
  # A pathological case where every bootstrap discrepancy could be smaller
  # than observed: mean(null_stats >= observed) would be exactly 0 under
  # the old formula. p must always be a multiple of 1/(B + 1), at least
  # 1/(B + 1).
  good <- ppe_dofv(n = 400L, df = 1, ncp = 9, seed = 91L)
  fit <- .ppeChiSquareMle(good, df = 1)

  p <- .ppeDiagnosticPValue(
    fit$retained, df = 1, estimate = fit$estimate, target = "ncp",
    bootstrapSamples = 10L, seed = 5L
  )

  expect_equal(p * 11, round(p * 11))
  expect_gte(p, 1 / 11)
})

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

test_that("a single call with two comparisons facets per-comparison, not collapsed", {
  # A genuine multi-comparison call, not two single-comparison calls glued
  # together afterward: fake_ppe_mixed_sse_object() is a three-model fixture
  # (simulation, alt1, alt2) built so a power comparison (simulation vs alt1)
  # and a Type-I comparison (alt2 vs simulation) are BOTH statistically valid
  # in the same run, letting facet_wrap()'s length(cmps) > 1 branch actually
  # execute rather than sit untested.
  skip_if_not_installed("patchwork")
  sse <- fake_ppe_mixed_sse_object(df = 1, powerNcp = 10, n = 80L, seed = 44L)
  power_cmp <- sseComparison("simulation", "alt1", df = 1, label = "power check")
  type1_cmp <- sseComparison("alt2", "simulation", df = 1, label = "type1 check")

  p <- plotSSEPpeDiagnostics(
    sse, comparisons = list(power_cmp, type1_cmp),
    bootstrapSamples = 20L, bootSeed = 5L
  )

  diag <- attr(p, "ppeDiagnostics")
  expect_equal(nrow(diag), 2L)
  expect_equal(diag$comparison, c("power check", "type1 check"))
  expect_equal(diag$parameter, c("ncp", "df"))
  # Each row's estimate is close to the true value that generated it (ncp=10,
  # df=1), confirming the two rows are not both diagnosing the same fit.
  expect_equal(diag$estimate[diag$parameter == "ncp"], 10, tolerance = 0.5)
  expect_equal(diag$estimate[diag$parameter == "df"], 1, tolerance = 0.5)

  # And the facets themselves are comparison-specific, not one curve shared
  # across both panels: each ECDF layer's points belong to exactly one facet.
  built <- ggplot2::ggplot_build(p$patches$plots[[1]])
  ecdf_layer <- built$data[[2]]
  panel_by_group <- table(PANEL = ecdf_layer$PANEL, group = ecdf_layer$group)
  expect_equal(unname(rowSums(panel_by_group > 0)), c(1L, 1L))
})

test_that("a criticalValue-only comparison is refused", {
  skip_if_not_installed("patchwork")
  sse <- fake_ppe_sse_object(df = 1, ncp = 10, n = 40L, seed = 45L)
  cmp <- sseComparison("simulation", "alt1", criticalValue = 5)

  expect_error(plotSSEPpeDiagnostics(sse, comparisons = cmp), "noncentral chi-square")
})
