skip_on_cran()

# ---- plotSSEParameterDraws() -----------------------------------------------

test_that("plotSSEParameterDraws returns a ggplot", {
  sse <- fake_draw_sse_object(mode = "joint", seed = 71L)
  p <- plotSSEParameterDraws(sse)

  expect_s3_class(p, "ggplot")
})

test_that("it facets one panel per drawn parameter", {
  # A real faceting check via ggplot2::ggplot_build(), not just "it runs" --
  # mirrors the technique used in test-plot-sse-ppe-diagnostics.R's
  # multi-comparison facet test.
  sse <- fake_draw_sse_object(mode = "joint", seed = 72L)
  p <- plotSSEParameterDraws(sse)

  built <- ggplot2::ggplot_build(p)
  n_panels <- nrow(unique(built$data[[1]]["PANEL"]))
  n_parameters <- length(unique(sse$initialValues$parameter))

  expect_equal(n_parameters, 2L)
  expect_equal(n_panels, n_parameters)
})

test_that("it also facets correctly for independent_iw draws", {
  sse <- fake_draw_sse_object(mode = "independent_iw", seed = 73L)
  p <- plotSSEParameterDraws(sse)

  built <- ggplot2::ggplot_build(p)
  n_panels <- nrow(unique(built$data[[1]]["PANEL"]))

  expect_equal(n_panels, 2L)
})

test_that("the parameterDrawSummary attribute matches a direct call", {
  sse <- fake_draw_sse_object(mode = "joint", seed = 74L)
  p <- plotSSEParameterDraws(sse)

  expect_equal(attr(p, "parameterDrawSummary"), parameterDrawSummary(sse))
})

test_that("a non-covariance-mode run errors via parameterDrawSummary()'s own message", {
  sse <- fake_ppe_sse_object(df = 1, ncp = 10, n = 40L, seed = 75L)

  expect_error(plotSSEParameterDraws(sse), "covariance")
})
