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
  # dispersion_ratio's own formula, restated -- would pass even if
  # realized_sd itself were wrong. The concrete check below is what actually
  # pins the value: the fixture's single eta IS the binding element, so its
  # realized dispersion should track its reported target SD closely (ratio
  # near 1), unlike a non-binding element in a larger block, which would come
  # out systematically above 1 (see the "independent_iw" module docs).
  expect_equal(row$dispersion_ratio, row$realized_sd / row$target_sd)
  expect_equal(row$dispersion_ratio, 1, tolerance = 0.25)
})

test_that("joint mode reports pairwise empirical dependence; independent_iw omits it", {
  # Code-quality review finding: the "jointDependence" attribute (empirical
  # correlation between jointly-drawn parameters, reported descriptively --
  # never claimed to reproduce fit$cov exactly) had zero test coverage.
  joint_sse <- fake_draw_sse_object(mode = "joint", omega = 0.6, se = 0.3,
                                    n = 300L, seed = 71L)
  s_joint <- parameterDrawSummary(joint_sse)
  dep <- attr(s_joint, "jointDependence")

  expect_true(is.data.frame(dep))
  expect_equal(nrow(dep), 1L)
  expect_setequal(
    c(dep$parameter1, dep$parameter2),
    c("tka", "omega(eta.ka,eta.ka)")
  )
  expect_true(is.finite(dep$correlation))
  expect_true(dep$correlation >= -1 && dep$correlation <= 1)

  # covarianceDraw = "independent_iw" draws THETA and OMEGA independently, so
  # there is no joint dependence structure to report at all -- the attribute
  # is absent, not an empty/NA table.
  iw_sse <- fake_draw_sse_object(mode = "independent_iw", omega = 0.6, se = 0.3,
                                 n = 300L, seed = 72L)
  s_iw <- parameterDrawSummary(iw_sse)
  expect_null(attr(s_iw, "jointDependence"))
})

test_that("parameterDrawSummary refuses a non-covariance-mode run", {
  sse <- fake_ppe_sse_object(n = 20L, seed = 66L)
  expect_error(parameterDrawSummary(sse), "covariance")
})

test_that("parameterDrawSummary refuses a non-nlmixr2SSE input", {
  expect_error(parameterDrawSummary(list()), "nlmixr2SSE")
})
