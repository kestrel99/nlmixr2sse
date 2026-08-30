skip_on_cran()

test_that("bias MCSE is the SD of replicate errors over sqrt(n)", {
  errors <- c(0.1, -0.2, 0.3, 0.05)
  expect_equal(.mcseFromErrors(errors), stats::sd(errors) / sqrt(4))
})

test_that("MCSE is nonnegative for a negative fixed truth", {
  s <- .parameterSummaryRow(estimates = c(-1.8, -2.2, -2.1), truth = -2)

  expect_gte(s$mcse_bias, 0)
  expect_gte(s$mcse_relative_bias, 0)
})

test_that("zero truths are excluded from relative bias", {
  s <- .parameterSummaryRow(estimates = c(0.1, 0.2), truth = 0)

  expect_true(is.na(s$mcse_relative_bias))
  expect_equal(s$n_effective_relative, 0L)
})

test_that("a single replicate yields NA intervals with its count retained", {
  s <- .parameterSummaryRow(estimates = 1.5, truth = 1)

  expect_true(is.na(s$ci_bias_lower))
  expect_equal(s$n_effective, 1L)
})

test_that(".parameterSummaryRow matches a hand-computed example", {
  # A fully worked example so a swapped formula (e.g. population SD, or the
  # relative branch dividing by something other than truth) would fail this
  # even though it might still pass the sign/edge-case tests above.
  estimates <- c(9, 11, 8, 12)
  truth <- 10
  errors <- estimates - truth
  rel_errors <- errors / truth

  s <- .parameterSummaryRow(estimates, truth)

  expect_equal(s$mcse_bias, stats::sd(errors) / sqrt(4))
  expect_equal(s$mcse_relative_bias, stats::sd(rel_errors) / sqrt(4))
  expect_equal(s$n_effective, 4L)
  expect_equal(s$n_effective_relative, 4L)

  z <- stats::qnorm(0.975)
  bias <- mean(errors)
  expect_equal(s$ci_bias_lower, bias - z * s$mcse_bias)
  expect_equal(s$ci_bias_upper, bias + z * s$mcse_bias)
})

test_that("rse is a superseded alias equal to mcse_relative_bias, not the old signed formula", {
  # Regression test for the historical bug: the old `rse` formula divided the
  # (nonnegative) SD of raw estimates by the truth directly, so a negative
  # fixed truth produced a NEGATIVE "standard error" -- nonsensical for a
  # quantity that is, by definition, a standard error. `.mcseFromErrors()`'s
  # `sd()/sqrt(n)` is nonnegative by construction regardless of truth's sign,
  # so both `rse` and `mcse_relative_bias` must come out nonnegative and
  # identical to each other here.
  raw_results <- data.frame(
    sample = c(1L, 2L, 3L),
    model_label = "ref_fit",
    role = "simulation",
    error_message = NA_character_,
    tka = c(-1.8, -2.2, -2.1),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  attr(raw_results, "rawResultsHeader") <- list(
    schema_version = 1L,
    columns = names(raw_results),
    base_cols = names(raw_results),
    theta_cols = "tka",
    omega_cols = character(0),
    sigma_cols = character(0),
    se_cols = character(0),
    parameter_cols = "tka"
  )
  initial_wide <- data.frame(
    sample = c(1L, 2L, 3L),
    tka = rep(-2, 3),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  spec_map <- list(ref_fit = list(
    thetaCols = "tka", omegaCols = character(0), sigmaCols = character(0)
  ))

  summary_tbl <- .computeParameterSummary(
    rawResults = raw_results, initialWide = initial_wide, specMap = spec_map
  )

  rse_row <- subset(summary_tbl, parameter == "tka" & statistic == "rse")
  mcse_row <- subset(
    summary_tbl, parameter == "tka" & statistic == "mcse_relative_bias"
  )

  expect_gte(rse_row$value, 0)
  expect_equal(rse_row$value, mcse_row$value)
})
