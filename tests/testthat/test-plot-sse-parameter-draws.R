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

test_that("it also facets correctly and places the reference line for independent_iw draws", {
  sse <- fake_draw_sse_object(mode = "independent_iw", omega = 0.6, seed = 73L)
  p <- plotSSEParameterDraws(sse)
  ds <- parameterDrawSummary(sse)

  built <- ggplot2::ggplot_build(p)
  n_panels <- nrow(unique(built$data[[1]]["PANEL"]))
  expect_equal(n_panels, 2L)

  panel_parameter <- built$layout$layout[, c("PANEL", "parameter")]
  vline_layer <- built$data[[2]]
  actual <- merge(vline_layer, panel_parameter, by = "PANEL")
  actual <- actual[order(actual$parameter), ]
  expected <- ds[order(ds$parameter), c("parameter", "target_mean")]

  expect_equal(actual$xintercept, expected$target_mean)
})

test_that("the reference line sits at each panel's actual target_mean", {
  # Code-quality review finding: prior tests confirmed the plot facets and
  # carries the right attribute, but never checked the geom_vline layer's
  # real xintercept against target_mean -- the one property that makes the
  # line correct rather than merely present. Mirrors the layer-data
  # inspection technique used in test-plot-sse-ppe-diagnostics.R.
  sse <- fake_draw_sse_object(mode = "joint", omega = 0.6, seed = 76L)
  p <- plotSSEParameterDraws(sse)
  ds <- parameterDrawSummary(sse)

  built <- ggplot2::ggplot_build(p)
  panel_parameter <- built$layout$layout[, c("PANEL", "parameter")]
  vline_layer <- built$data[[2]]

  actual <- merge(vline_layer, panel_parameter, by = "PANEL")
  actual <- actual[order(actual$parameter), ]
  expected <- ds[order(ds$parameter), c("parameter", "target_mean")]

  expect_equal(actual$xintercept, expected$target_mean)
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
