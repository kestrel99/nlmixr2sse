# Visual counterpart to parameterDrawSummary() (R/sse-parameter-diagnostics.R).
#
# parameterDrawSummary() reports, numerically, how well realized replicate
# parameter draws from a `parameterSource = "covariance"` run match the
# targets that seeded them (the fitted value and reported SE from
# `x$runInfo$parameterSourceInfo$targets`) -- including the two documented
# draw-mode artefacts described in that file's own header (log-Cholesky mean
# drift for `covarianceDraw = "joint"`, non-binding-element over-dispersion
# for `covarianceDraw = "independent_iw"`). This file adds the visual
# equivalent: one histogram panel per drawn parameter, so those artefacts are
# visible as a shape, not just a number.

#' Plot realized covariance-mode parameter draws against their targets
#'
#' For a completed `runSSE()` run with `parameterSource = "covariance"`, this
#' plots -- for every drawn generating parameter (THETA or OMEGA) -- a
#' histogram of the ACTUALLY DRAWN replicate values (`x$initialValues`), with
#' a vertical reference line at that parameter's `target_mean` (the fitted
#' value the draw was seeded from, from `parameterDrawSummary(x)`). It is the
#' visual counterpart to `parameterDrawSummary()`, which reports the same
#' comparison numerically -- see that function's documentation for what the
#' two known draw-mode artefacts (`covarianceDraw = "joint"` raw-scale mean
#' drift and `covarianceDraw = "independent_iw"` non-binding-element
#' over-dispersion) look like in the histogram: a joint-mode histogram
#' visibly shifted right of its reference line, or an independent_iw-mode
#' histogram visibly wider than its target SE would suggest.
#'
#' All validation is delegated to `parameterDrawSummary(x)`, called first: it
#' asserts `x` is an `nlmixr2SSE` object and aborts if
#' `x$runInfo$parameterSource` is not `"covariance"`. This function does not
#' duplicate that check or its error message.
#'
#' @param x An `nlmixr2SSE` object produced by a `parameterSource =
#'   "covariance"` `runSSE()` call.
#' @param ... Reserved for future extension.
#' @return A `ggplot` object, one facet panel per drawn parameter
#'   (`ggplot2::facet_wrap(~parameter, scales = "free")`), with
#'   `parameterDrawSummary(x)`'s return value attached as a
#'   `"parameterDrawSummary"` attribute -- so the numeric diagnostic (target
#'   means/SDs, `mean_drift_flag`, `dispersion_ratio`, and the rest of that
#'   table) is always available without reading pixels off the plot.
#' @export
plotSSEParameterDraws <- function(x, ...) {
  .assertNamespace("ggplot2", "construct SSE plots")
  draw_summary <- parameterDrawSummary(x)

  long <- x$initialValues
  targets <- draw_summary[, c("parameter", "target_mean")]

  p <- ggplot2::ggplot(long, ggplot2::aes(x = value)) +
    # 20 bins, matching plotSSEOfvDistribution()'s histogram convention
    # (R/plot-sse.R) -- an explicit choice rather than ggplot2's bins = 30
    # default, which is more likely to fragment a small replicate count into
    # a sparse, hard-to-read histogram.
    ggplot2::geom_histogram(bins = 20L) +
    ggplot2::geom_vline(
      data = targets,
      ggplot2::aes(xintercept = target_mean),
      colour = "firebrick", linetype = 2
    ) +
    ggplot2::facet_wrap(~parameter, scales = "free") +
    ggplot2::labs(
      x = "Realized draw value", y = "Count",
      title = "Realized parameter draws vs. target mean",
      subtitle = paste(
        "Dashed line: target_mean",
        "(the fitted value the draw was seeded from)"
      )
    ) +
    ggplot2::theme_minimal()

  attr(p, "parameterDrawSummary") <- draw_summary
  p
}
