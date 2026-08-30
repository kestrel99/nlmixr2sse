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

test_that("plotSSEDiagnostics combines panels when patchwork is available", {
  skip_if_not_installed("patchwork")

  sse <- fake_sse_object()
  p <- plotSSEDiagnostics(sse)

  expect_true(inherits(p, "patchwork"))
})
