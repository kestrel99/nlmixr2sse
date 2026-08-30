skip_on_cran()

test_that("plotSSEParameterBias returns a ggplot backed by parameter summary data", {
  sse <- fake_sse_object()

  p <- plotSSEParameterBias(
    sse,
    statistic = "relative_bias"
  )

  expect_s3_class(p, "ggplot")
  expect_equal(
    p$data,
    .parameterBiasPlotData(sse, statistic = "relative_bias")
  )
})

test_that("plotSSEParameterEstimates returns a ggplot backed by replicate-level data", {
  sse <- fake_sse_object()

  p <- plotSSEParameterEstimates(sse)

  expect_s3_class(p, "ggplot")
  expect_equal(
    p$data,
    .parameterEstimatePlotData(sse)
  )
})

test_that("plotSSEOfvDistribution returns a ggplot backed by delta OFV data", {
  sse <- fake_sse_object()

  p <- plotSSEOfvDistribution(sse)

  expect_s3_class(p, "ggplot")
  expect_equal(
    p$data,
    .ofvDeltaPlotData(sse)
  )
})

test_that("plotSSEPower returns a ggplot backed by OFV summary data", {
  sse <- fake_sse_object()

  p <- plotSSEPower(sse, direction = "power")

  expect_s3_class(p, "ggplot")
  expect_equal(
    p$data,
    .powerPlotData(sse, direction = "power")
  )
})

test_that("plotSSEPpePower returns a ggplot backed by PPE curve data", {
  sse <- fake_sse_object()

  # fake_sse_object() has no explicit sseComparison() in runInfo$comparisons,
  # so the distribution_mle default infers df from parameter counts and
  # (deliberately, since Task 5 restored ppe = TRUE) warns about it. That
  # warning is exercised on its own below; here it is incidental to what
  # this test checks, so it is suppressed rather than asserted twice.
  p <- suppressWarnings(
    plotSSEPpePower(sse, thresholds = 3.84, studySizes = c(6L, 12L, 18L))
  )

  expect_s3_class(p, "ggplot")
  expect_equal(
    p$data,
    suppressWarnings(.ppePowerPlotData(
      sse,
      thresholds = 3.84,
      studySizes = c(6L, 12L, 18L)
    ))
  )
})

test_that("plotSSEPpePower warns about inferred df exactly once, not zero or twice", {
  # Regression test: plotSSEPpePower() resolves comparisons itself (to split
  # by mode) and then calls .ppePowerPlotData(), which resolves again inside
  # .ppePowerPlotDataMle(). A version that resolved-and-suppressed for the
  # mode split, then relied on the second "real" pass to raise the
  # inferred-df warning, silently warned ZERO times: handing an
  # already-resolved (non-NULL) comparisons list to the second pass means it
  # never re-enters .resolveComparisons()'s comparisons-is-NULL branch, the
  # only place that warning is ever raised. Calling .ppePowerPlotData()
  # directly (as the test above does) resolves only once and does not
  # exercise this two-pass path at all -- only plotSSEPpePower() itself can
  # catch a regression here.
  sse <- fake_sse_object()

  warnings <- testthat::capture_warnings(
    plotSSEPpePower(sse, thresholds = 3.84, studySizes = c(6L, 12L, 18L))
  )

  expect_equal(sum(grepl("Degrees of freedom inferred", warnings)), 1L)
})

test_that("plot.nlmixr2SSE dispatches to named plot helpers", {
  sse <- fake_sse_object()

  expect_s3_class(plot(sse, type = "parameter_bias"), "ggplot")
  expect_s3_class(plot(sse, type = "parameter_estimates"), "ggplot")
  expect_s3_class(plot(sse, type = "ofv_distribution"), "ggplot")
  expect_s3_class(plot(sse, type = "power"), "ggplot")
  # See the comment above: fake_sse_object() has no explicit comparison, so
  # the distribution_mle default warns about its inferred df.
  expect_s3_class(suppressWarnings(plot(sse, type = "ppe_power")), "ggplot")
})

test_that("a power comparison's ribbon is a real bootstrap interval, not still NA", {
  sse <- fake_ppe_sse_object(df = 1, ncp = 10, n = 100L, seed = 61L)
  cmp <- sseComparison("simulation", "alt1", df = 1)

  data <- .ppePowerPlotData(
    sse, thresholds = 3.84, studySizes = c(6L, 12L, 18L),
    comparisons = cmp, nonpositivePolicy = "drop",
    bootstrapSamples = 100L, bootSeed = 7L
  )

  expect_true(all(is.finite(data$power_lower)))
  expect_true(all(is.finite(data$power_upper)))
  # The interval must bracket the point estimate for at least one row, or a
  # "real" ribbon that happens to sit on top of the curve would pass this
  # check for the wrong reason.
  expect_true(any(data$power_lower < data$power & data$power < data$power_upper))
})

test_that("a type1 comparison renders a point-range, not a curve", {
  sse <- fake_ppe_type1_sse_object(df = 1, n = 200L, seed = 71L)
  p <- plotSSEPpePower(
    sse, comparisons = sseComparison("alt1", "simulation", df = 1),
    bootstrapSamples = 50L, bootSeed = 3L
  )

  geoms <- vapply(p$layers, function(l) class(l$geom)[[1L]], character(1))
  expect_true("GeomPointrange" %in% geoms)
  expect_false("GeomLine" %in% geoms)
})

test_that("mixing power and type1 comparisons warns rather than plotting both", {
  sse <- fake_ppe_sse_object(df = 1, ncp = 8, n = 100L, seed = 72L)
  expect_warning(
    plotSSEPpePower(sse, comparisons = list(
      sseComparison("simulation", "alt1", df = 1),
      sseComparison("alt1", "simulation", df = 1, label = "type1 check")
    ), bootstrapSamples = 20L, bootSeed = 3L),
    "omitted"
  )
})

test_that("mixing power and type1 comparisons still renders the power curve", {
  sse <- fake_ppe_sse_object(df = 1, ncp = 8, n = 100L, seed = 73L)
  p <- suppressWarnings(plotSSEPpePower(sse, comparisons = list(
    sseComparison("simulation", "alt1", df = 1),
    sseComparison("alt1", "simulation", df = 1, label = "type1 check")
  ), bootstrapSamples = 20L, bootSeed = 3L))

  geoms <- vapply(p$layers, function(l) class(l$geom)[[1L]], character(1))
  expect_true("GeomLine" %in% geoms)
  expect_false("GeomPointrange" %in% geoms)
})

test_that("plotSSEDiagnostics combines panels when patchwork is available", {
  skip_if_not_installed("patchwork")

  sse <- fake_sse_object()
  p <- plotSSEDiagnostics(sse)

  expect_true(inherits(p, "patchwork"))
})
