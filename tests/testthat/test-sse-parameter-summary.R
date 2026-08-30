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

test_that("the long-format table applies the percentage scaling exactly once, to exactly the relative fields", {
  # .parameterSummaryRow() itself returns fractional (0-1) values; the
  # integration point in .singleParameterSummary() multiplies only the
  # relative-scale fields by 100 to match relative_bias/rse's existing
  # percentage convention, leaving the absolute-scale fields (mcse_bias,
  # ci_bias_*) untouched. A systematic scaling bug here (an extra 100x, a
  # missing one, or scaling the wrong fields) would not be "obviously wrong"
  # -- the self-consistency check above (rse == mcse_relative_bias) cannot
  # catch it, since both sides share the same scaling code path. Pin
  # everything to hand-computed values instead.
  estimates <- c(9, 12, 8, 13)
  truth <- 10
  errors <- estimates - truth
  rel_errors <- errors / truth
  z <- stats::qnorm(0.975)

  raw_results <- data.frame(
    sample = 1:4,
    model_label = "ref_fit",
    role = "simulation",
    error_message = NA_character_,
    tka = estimates,
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
    sample = 1:4,
    tka = rep(truth, 4),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  spec_map <- list(ref_fit = list(
    thetaCols = "tka", omegaCols = character(0), sigmaCols = character(0)
  ))

  summary_tbl <- .computeParameterSummary(
    rawResults = raw_results, initialWide = initial_wide, specMap = spec_map
  )
  value_of <- function(stat) {
    subset(summary_tbl, parameter == "tka" & statistic == stat)$value
  }

  expect_equal(value_of("mcse_bias"), stats::sd(errors) / sqrt(4))
  expect_equal(value_of("ci_bias_lower"), mean(errors) - z * stats::sd(errors) / sqrt(4))
  expect_equal(value_of("ci_bias_upper"), mean(errors) + z * stats::sd(errors) / sqrt(4))

  expect_equal(value_of("mcse_relative_bias"), 100 * stats::sd(rel_errors) / sqrt(4))
  expect_equal(
    value_of("ci_relative_bias_lower"),
    100 * (mean(rel_errors) - z * stats::sd(rel_errors) / sqrt(4))
  )
  expect_equal(
    value_of("ci_relative_bias_upper"),
    100 * (mean(rel_errors) + z * stats::sd(rel_errors) / sqrt(4))
  )
})
