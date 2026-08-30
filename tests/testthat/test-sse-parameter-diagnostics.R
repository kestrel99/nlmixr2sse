skip_on_cran()

test_that("joint mode exposes raw-scale OMEGA mean drift", {
  sse <- fake_draw_sse_object(mode = "joint", omega = 0.6, se = 0.3,
                              n = 500L, seed = 61L)
  s <- parameterDrawSummary(sse)
  row <- s[s$parameter == "omega(eta.ka,eta.ka)", ]

  # The log-Cholesky back-transformation inflates the raw mean by
  # exp(SE^2 / (2 * Omega^2)); at SE/Omega = 0.5 that is about 13%.
  expect_gt(row$realized_mean / row$target_mean, 1.05)
  expect_true(row$mean_drift_flag)
})

test_that("independent_iw reports binding degrees of freedom per block", {
  sse <- fake_draw_sse_object(mode = "independent_iw", omega = 0.6, se = 0.3,
                              n = 200L, seed = 62L)
  s <- parameterDrawSummary(sse)

  expect_true(all(s$binding_nu[!is.na(s$binding_nu)] > 0))
})

test_that("every drawn OMEGA block is positive definite", {
  sse <- fake_draw_sse_object(mode = "joint", omega = 0.6, se = 0.3,
                              n = 200L, seed = 63L)
  s <- parameterDrawSummary(sse)

  expect_equal(sum(s$n_not_positive_definite, na.rm = TRUE), 0L)
})

test_that("out-of-bound THETA draws are counted, never resampled", {
  sse <- fake_draw_sse_object(mode = "joint", thetaLower = 0.5,
                              n = 200L, seed = 64L)
  s <- parameterDrawSummary(sse)

  expect_gt(sum(s$n_out_of_domain, na.rm = TRUE), 0L)
  expect_equal(nrow(s), length(unique(s$parameter)))

  # The underlying draws themselves are never truncated or resampled -- the
  # out-of-domain count is diagnostic only.
  raw_tka <- sse$initialValues$value[sse$initialValues$parameter == "tka"]
  expect_true(any(raw_tka < 0.5))
})

test_that("realized dispersion summary statistics are reported per parameter", {
  sse <- fake_draw_sse_object(mode = "independent_iw", omega = 0.6, se = 0.3,
                              n = 300L, seed = 65L)
  s <- parameterDrawSummary(sse)
  row <- s[s$parameter == "omega(eta.ka,eta.ka)", ]

  expect_true(all(c(
    "realized_mean", "realized_sd", "realized_median",
    "realized_q025", "realized_q975", "realized_min", "realized_max", "n",
    "dispersion_ratio"
  ) %in% names(s)))
  expect_equal(row$n, 300L)
  expect_true(row$realized_min <= row$realized_median)
  expect_true(row$realized_median <= row$realized_max)
  # The binding element's realized SD should track its reported target SD
  # reasonably closely (this fixture's single eta IS the binding element).
  expect_equal(row$dispersion_ratio, row$realized_sd / row$target_sd)
})

test_that("parameterDrawSummary refuses a non-covariance-mode run", {
  sse <- fake_ppe_sse_object(n = 20L, seed = 66L)
  expect_error(parameterDrawSummary(sse), "covariance")
})

test_that("parameterDrawSummary refuses a non-nlmixr2SSE input", {
  expect_error(parameterDrawSummary(list()), "nlmixr2SSE")
})
