# These tests pin down TODAY's statistical behaviour of the SSE outputs so
# that later tasks -- which deliberately change model comparisons, paired
# denominators, and the PPE estimator -- make that change visibly rather than
# silently. They are expected to PASS now and to FAIL only when a later task
# intentionally alters the behaviour they describe.

skip_on_cran()

test_that("the paired fixture reproduces the supplied OFVs exactly", {
  sse <- fake_paired_sse_object(full_ofv = c(100, 101), reduced_ofv = c(105, 108))
  expect_equal(-.ofvDeltaPlotData(sse)$delta_ofv, c(5, 7))
})

test_that("the PPE fixture reproduces its draws as test statistics", {
  d <- ppe_dofv(n = 20L, df = 1, ncp = 8, seed = 5L)
  sse <- fake_ppe_sse_object(df = 1, ncp = 8, n = 20L, seed = 5L)
  expect_equal(-.ofvDeltaPlotData(sse)$delta_ofv, d)
})

test_that("ppe_dofv leaves the caller's RNG untouched", {
  set.seed(1L)
  before <- .Random.seed
  ppe_dofv(n = 10L, seed = 2L)
  expect_equal(.Random.seed, before)
})

test_that("current delta sign convention is reference minus alternative", {
  sse <- fake_sse_object()
  delta <- .ofvDeltaPlotData(sse)

  # simulation objf 100/101, alternative objf 103/105
  expect_equal(delta$delta_ofv, c(-3, -4))
  # and the PPE test statistic is the negation
  expect_equal(-delta$delta_ofv, c(3, 4))
})

test_that("current PPE inverts the exceedance probability per threshold", {
  sse <- fake_sse_object()

  low <- .ppePowerPlotData(sse, thresholds = 1, studySizes = 12L)
  high <- .ppePowerPlotData(sse, thresholds = 3.5, studySizes = 12L)

  # a separate ncp per threshold is exactly the behaviour being replaced
  expect_false(identical(low$power, high$power))
  expect_equal(unique(low$degrees_freedom), 1L)
})

test_that("modelDegreesFreedom falls back to 1 for an unknown label", {
  sse <- fake_sse_object()
  expect_equal(.modelDegreesFreedom(sse, "no_such_model"), 1L)
})
