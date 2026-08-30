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

test_that("ppe_dofv restores the absence of .Random.seed", {
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    rm(".Random.seed", envir = globalenv())
  }
  ppe_dofv(n = 5L, seed = 3L)
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
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
  test_stat <- c(3, 4)  # -delta_ofv for the two-sample fixture

  # Task 5 made "distribution_mle" (a single whole-distribution noncentrality)
  # the default; this test pins the older per-threshold estimator, still
  # available as method = "exceedance", so it must ask for it explicitly.
  low <- .ppePowerPlotData(sse, thresholds = 1, studySizes = 12L, method = "exceedance")
  high <- .ppePowerPlotData(sse, thresholds = 3.5, studySizes = 12L, method = "exceedance")

  # At the base study size the scaling factor is 1, so each threshold's ncp is
  # solved to reproduce that threshold's own clipped exceedance rate exactly.
  # A single ncp fitted to the whole distribution cannot reproduce two
  # different rates at once, so these equalities are precisely what breaks when
  # the estimator changes in Task 5. Asserting only that the two powers differ
  # would be useless: power depends on the threshold under any implementation.
  expect_equal(low$power, 100 * .clipProbability(mean(test_stat > 1), 2))
  expect_equal(high$power, 100 * .clipProbability(mean(test_stat > 3.5), 2))
  expect_equal(unique(low$degrees_freedom), 1L)
})

test_that("modelDegreesFreedom falls back to 1 for an unknown label", {
  sse <- fake_sse_object()
  expect_equal(.modelDegreesFreedom(sse, "no_such_model"), 1L)
})

test_that("fake_paired_sse_object requires at least one finite OFV", {
  expect_error(
    fake_paired_sse_object(full_ofv = numeric(0), reduced_ofv = numeric(0)),
    "at least one finite OFV"
  )
  expect_error(
    fake_paired_sse_object(full_ofv = NA_real_, reduced_ofv = NA_real_),
    "at least one finite OFV"
  )
})

test_that("fake_ppe_sse_object carries the real seed into runInfo", {
  sse <- fake_ppe_sse_object(n = 5L, seed = 7L)
  expect_equal(sse$runInfo$seed, 7L)
})

test_that("fake_ppe_sse_object rejects nNonPositive greater than n", {
  expect_error(fake_ppe_sse_object(n = 5L, nNonPositive = 9L))
})
