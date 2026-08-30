test_that("the denominator is paired evaluable replicates, not any single model", {
  sse <- fake_paired_sse_object(
    full_ofv    = c(100, 101, 102, NA,  104),
    reduced_ofv = c(110, 111, NA,  113, 114)
  )
  cmp <- sseComparison("simulation", "alt1", df = 1)

  s <- comparisonSummary(sse, comparisons = cmp)

  expect_equal(s$n_attempted, 5L)
  expect_equal(s$n_full_evaluable, 4L)
  expect_equal(s$n_reduced_evaluable, 4L)
  expect_equal(s$n_paired_evaluable, 3L)   # samples 1, 2, 5
  expect_equal(s$n_excluded, 2L)
})

test_that("the empirical interval is reproducible from the reported counts", {
  sse <- fake_paired_sse_object(
    full_ofv    = rep(100, 20),
    reduced_ofv = 100 + c(rep(10, 12), rep(0.5, 8))
  )
  cmp <- sseComparison("simulation", "alt1", df = 1)

  s <- comparisonSummary(sse, comparisons = cmp)

  expect_equal(s$n_exceeding, 12L)
  expect_equal(s$probability, 12 / 20)
  expect_equal(
    c(s$ci_lower, s$ci_upper),
    c(stats::qbeta(0.025, 12, 9), stats::qbeta(0.975, 13, 8))
  )
})

test_that("the test statistic is reduced minus full regardless of which was simulated", {
  sse <- fake_paired_sse_object(full_ofv = c(100, 100), reduced_ofv = c(105, 108))

  power_cmp <- sseComparison("simulation", "alt1", df = 1)
  expect_equal(.comparisonTestStatistic(sse, .resolveComparison(sse, power_cmp)), c(5, 8))
})

# The plan's tests never actually exercise a Type-I comparison (simulation
# model is the reduced one) through comparisonSummary(). That is the whole
# point of computing the statistic from the comparison rather than a sign
# convention -- .ofvDeltaPlotData() gets this backwards for exactly this case
# -- so it needs direct coverage, not just an inference from the power case.
test_that("a Type-I comparison yields mode 'type1' and a positive test statistic", {
  # fake_paired_sse_object()'s `full_ofv`/`reduced_ofv` args just name the
  # OFV vectors assigned to the simulation model and "alt1" respectively --
  # they are not the "full"/"reduced" of the comparison under test. Here the
  # simulation model is the *reduced* member of the comparison (mode
  # "type1"), so nesting requires its OFV to be >= the more flexible alt1
  # model's OFV: simulation (full_ofv) = 105, 108; alt1 (reduced_ofv) = 100, 100.
  sse <- fake_paired_sse_object(full_ofv = c(105, 108), reduced_ofv = c(100, 100))
  cmp <- sseComparison("alt1", "simulation", df = 1)

  s <- comparisonSummary(sse, comparisons = cmp)

  expect_equal(s$mode, "type1")
  stat <- .comparisonTestStatistic(sse, .resolveComparison(sse, cmp))
  expect_true(all(stat > 0))
  expect_equal(stat, c(5, 8))
})

test_that("the minPairedFraction warning fires below the threshold and not above it", {
  # 3 of 5 attempted samples are paired evaluable: 60%.
  sse <- fake_paired_sse_object(
    full_ofv    = c(100, 101, 102, NA,  104),
    reduced_ofv = c(110, 111, NA,  113, 114)
  )
  cmp <- sseComparison("simulation", "alt1", df = 1)

  expect_warning(
    comparisonSummary(sse, comparisons = cmp, minPairedFraction = 0.7),
    "60%"
  )
  expect_no_warning(
    comparisonSummary(sse, comparisons = cmp, minPairedFraction = 0.5)
  )
})

test_that("a comparison with zero paired evaluable replicates aborts, naming the comparison", {
  sse <- fake_paired_sse_object(
    full_ofv    = c(100, NA),
    reduced_ofv = c(NA, 110)
  )
  cmp <- sseComparison("simulation", "alt1", df = 1, label = "no_overlap")

  expect_error(
    comparisonSummary(sse, comparisons = cmp),
    "no_overlap"
  )
})

test_that("multiple comparisons produce one row each, in the order supplied", {
  sse <- fake_paired_sse_object(
    full_ofv    = c(100, 100, 100),
    reduced_ofv = c(105, 108, 110)
  )
  cmp1 <- sseComparison("simulation", "alt1", df = 1, label = "first")
  cmp2 <- sseComparison("alt1", "simulation", df = 1, label = "second")

  s <- comparisonSummary(sse, comparisons = list(cmp1, cmp2))

  expect_equal(nrow(s), 2L)
  expect_equal(s$comparison, c("first", "second"))
})

test_that("mcse_probability equals sqrt(p(1-p)/n) and is 0 at p = 0 or p = 1", {
  sse_mixed <- fake_paired_sse_object(
    full_ofv    = rep(100, 20),
    reduced_ofv = 100 + c(rep(10, 12), rep(0.5, 8))
  )
  s_mixed <- comparisonSummary(sse_mixed, comparisons = sseComparison("simulation", "alt1", df = 1))
  p <- 12 / 20
  expect_equal(s_mixed$mcse_probability, sqrt(p * (1 - p) / 20))

  # Every replicate exceeds the critical value: p = 1.
  sse_all <- fake_paired_sse_object(
    full_ofv    = rep(100, 10),
    reduced_ofv = rep(200, 10)
  )
  s_all <- comparisonSummary(sse_all, comparisons = sseComparison("simulation", "alt1", df = 1))
  expect_equal(s_all$probability, 1)
  expect_equal(s_all$mcse_probability, 0)

  # No replicate exceeds the critical value: p = 0.
  sse_none <- fake_paired_sse_object(
    full_ofv    = rep(100, 10),
    reduced_ofv = rep(100.001, 10)
  )
  s_none <- comparisonSummary(sse_none, comparisons = sseComparison("simulation", "alt1", df = 1))
  expect_equal(s_none$probability, 0)
  expect_equal(s_none$mcse_probability, 0)
})
