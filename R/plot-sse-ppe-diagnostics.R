# Distribution adequacy diagnostics for the PPE noncentral chi-square fit.
#
# ppeSummary()/plotSSEPpePower() (R/sse-ppe.R, R/plot-sse.R) fit a noncentral
# chi-square to the retained likelihood-ratio test statistics and report a
# point estimate (and bootstrap interval) from that fit. Nothing in that
# pipeline checks whether the noncentral chi-square is actually a good
# description of the data -- this file adds that check: a Cramer-von Mises
# discrepancy between the empirical CDF of the retained statistics and the
# FITTED chi-square CDF, a parametric-bootstrap p-value for it, and a visual
# ECDF/QQ diagnostic. This is evidence about approximation adequacy, not a
# pass/fail gate -- see plotSSEPpeDiagnostics()'s docs.

#' Cramer-von Mises discrepancy between retained statistics and a fitted chi-square
#'
#' `sum((F_n(x_i) - u_i)^2)`-style discrepancy computed via the standard
#' probability-integral-transform form: if `x` truly came from the fitted
#' distribution, `u = pchisq(x, df, ncp)` is (approximately) Uniform(0, 1),
#' and this is the ordinary one-sample Cramer-von Mises statistic for that
#' transformed sample against the Uniform(0, 1) reference.
#'
#' @param x Sorted or unsorted retained (positive) test statistics.
#' @param df,ncp The fitted distribution's degrees of freedom and
#'   noncentrality: whichever member is NOT the one being diagnosed is the
#'   comparison's known/fixed value, the other is the MLE estimate. See
#'   `target`.
#' @param target Which of `df`/`ncp` is the value being diagnosed --
#'   `"ncp"` for a power comparison (df fixed, ncp estimated), `"df"` for a
#'   Type-I comparison (ncp fixed at 0, df estimated). Mirrors
#'   `.ppeParametricBootstrap()`'s `target` argument. The discrepancy formula
#'   itself does not depend on `target` -- `df`/`ncp` are just whichever
#'   values describe the fitted distribution being checked -- but callers
#'   (`.ppeDiagnosticPValue()`, `plotSSEPpeDiagnostics()`) need it to know
#'   which parameter to re-estimate when redrawing under the null.
#' @return A single numeric discrepancy statistic (larger = worse fit).
#' @noRd
.ppeCramerVonMises <- function(x, df, ncp, target = c("ncp", "df")) {
  target <- match.arg(target)
  x <- sort(x)
  n <- length(x)
  u <- stats::pchisq(x, df = df, ncp = ncp)
  1 / (12 * n) + sum((u - (2 * seq_len(n) - 1) / (2 * n))^2)
}

#' Parametric-bootstrap p-value for the Cramer-von Mises discrepancy
#'
#' Draws `bootstrapSamples` synthetic samples from the FITTED noncentral
#' chi-square, refits `.ppeChiSquareMle()` to each (mirroring
#' `.ppeParametricBootstrap()`'s `target`-branched draw/refit exactly, so the
#' null distribution respects the same fixed/estimated roles the point
#' estimate itself was fit under), recomputes the Cramer-von Mises statistic
#' for each refit, and returns the plus-one-corrected proportion of the null
#' (bootstrap) discrepancies at least as large as the one observed --
#' `(1 + #{T_b >= T_obs}) / (B + 1)`, the standard Monte Carlo p-value
#' construction, which is never exactly zero even when no bootstrap
#' discrepancy reaches the observed one (the naive `mean(T_b >= T_obs)` can
#' report an impossible exact zero with a finite bootstrap sample). A small
#' p-value means the fitted distribution is a poor description of the real
#' data. A refit that errors is dropped from the null distribution rather
#' than aborting the whole bootstrap (mirrors `.ppeParametricBootstrap()`'s
#' `n_failed` handling).
#'
#' @param x Retained (positive) test statistics the original fit used.
#' @param df Degrees of freedom fixed by the comparison. Used only when
#'   `target = "ncp"` (ignored, but still a required argument, when
#'   `target = "df"`, where `estimate` itself is the df being diagnosed).
#' @param estimate The fitted point estimate from `.ppeChiSquareMle()`.
#' @param target `"ncp"` (power comparison) or `"df"` (Type-I comparison).
#' @param bootstrapSamples Number of bootstrap replicates.
#' @param seed Passed to `.withPpeSeed()`.
#' @return A single p-value in `(0, 1]`, computed with the standard
#'   Monte Carlo plus-one correction
#'   `(1 + sum(null_stats >= observed)) / (length(null_stats) + 1)` so it is
#'   never exactly zero regardless of bootstrap sample size, or `NA_real_`
#'   if every bootstrap refit failed.
#' @noRd
.ppeDiagnosticPValue <- function(x, df, estimate, target = c("ncp", "df"),
                                 bootstrapSamples = 1000L, seed = NULL) {
  target <- match.arg(target)
  ncp <- if (identical(target, "ncp")) estimate else 0
  fixed_df <- if (identical(target, "ncp")) df else estimate
  observed <- .ppeCramerVonMises(x, df = fixed_df, ncp = ncp, target = target)
  n <- length(x)
  null_stats <- .withPpeSeed(seed, {
    vapply(seq_len(bootstrapSamples), function(i) {
      draw <- if (identical(target, "ncp")) {
        stats::rchisq(n, df = df, ncp = estimate)
      } else {
        stats::rchisq(n, df = estimate, ncp = 0)
      }
      refit <- tryCatch(
        if (identical(target, "ncp")) {
          .ppeChiSquareMle(draw, df = df)$estimate
        } else {
          .ppeChiSquareMle(draw, ncp = 0)$estimate
        },
        error = function(e) NA_real_
      )
      if (is.na(refit)) {
        return(NA_real_)
      }
      refit_ncp <- if (identical(target, "ncp")) refit else 0
      refit_df <- if (identical(target, "ncp")) df else refit
      # draw[draw > 0]: refit statistics are drawn from rchisq(), which is
      # strictly positive with probability 1, so this mirrors (but does not
      # rely on) .ppeChiSquareMle()'s own retained-values truncation.
      .ppeCramerVonMises(draw[draw > 0], df = refit_df, ncp = refit_ncp, target = target)
    }, numeric(1))
  })
  null_stats <- null_stats[is.finite(null_stats)]
  if (length(null_stats) == 0L) {
    return(NA_real_)
  }
  (1 + sum(null_stats >= observed)) / (length(null_stats) + 1)
}

#' Pointwise bootstrap envelope for the ECDF panel
#'
#' Simulates `envelopeSamples` synthetic samples of size `n` from the FITTED
#' distribution held FIXED at its estimated parameter (unlike
#' `.ppeDiagnosticPValue()`'s null-distribution bootstrap, this does not
#' refit each draw -- the envelope answers "how much would the ECDF of a
#' correctly-sized sample from this exact fitted model wobble", which is a
#' question about sampling variability alone, not about re-estimation
#' uncertainty), computes each replicate's ECDF at the points in `grid`, and
#' takes pointwise `conf.level` quantiles across replicates.
#'
#' This is a POINTWISE envelope, not a simultaneous one: at any single grid
#' point it covers the stated fraction of replicate ECDF values, but the
#' chance that the real ECDF strays outside it at at least one of many grid
#' points is higher than `1 - conf.level`. `plotSSEPpeDiagnostics()`'s
#' caption says so explicitly.
#'
#' @param n Sample size to simulate per replicate (the number of retained
#'   statistics in the real data).
#' @param df,ncp The fitted distribution's parameters (fixed for every draw).
#' @param grid Numeric vector of x-values at which to evaluate each
#'   replicate's ECDF.
#' @param envelopeSamples Number of replicate samples. `0` returns `NA` bounds.
#' @param conf.level Confidence level for the pointwise envelope.
#' @param seed Passed to `.withPpeSeed()`.
#' @return A list with `lower` and `upper`, numeric vectors the same length
#'   as `grid`.
#' @noRd
.ppeEnvelopeBounds <- function(n, df, ncp, grid, envelopeSamples, conf.level, seed) {
  if (envelopeSamples <= 0L) {
    na <- rep(NA_real_, length(grid))
    return(list(lower = na, upper = na))
  }
  a <- 1 - conf.level
  ecdf_matrix <- .withPpeSeed(seed, {
    vapply(seq_len(envelopeSamples), function(i) {
      sim <- stats::rchisq(n, df = df, ncp = ncp)
      stats::ecdf(sim)(grid)
    }, numeric(length(grid)))
  })
  list(
    lower = apply(ecdf_matrix, 1L, stats::quantile, probs = a / 2, names = FALSE),
    upper = apply(ecdf_matrix, 1L, stats::quantile, probs = 1 - a / 2, names = FALSE)
  )
}

#' Compute the diagnostic table and panel data for one resolved comparison
#'
#' @param x An `nlmixr2SSE` object.
#' @param comparison A resolved `sseComparison`.
#' @param nonpositive,conf.level,bootstrapSamples,bootSeed As in
#'   `plotSSEPpeDiagnostics()`.
#' @param gridLength Number of x-values spanning the ECDF envelope grid.
#' @return A list with `row` (one-row `data.frame`), `ecdf_data`,
#'   `envelope_data`, and `qq_data` (each a `data.frame` tagged with
#'   `comparison`).
#' @noRd
.ppeDiagnosticsForComparison <- function(x, comparison, nonpositive, conf.level,
                                         bootstrapSamples, bootSeed,
                                         gridLength = 200L) {
  # Same target split .ppeFit() already uses: a power comparison's simulation
  # model is the full one, so ncp is estimated with df held at the
  # comparison's known value; a Type-I comparison's simulation model is the
  # reduced one, so df is estimated with ncp held at 0.
  target <- if (identical(comparison$mode, "type1")) "df" else "ncp"
  stat <- .comparisonTestStatistic(x, comparison)
  fit <- if (identical(target, "ncp")) {
    .ppeChiSquareMle(stat, df = comparison$df)
  } else {
    .ppeChiSquareMle(stat, ncp = 0)
  }
  .ppeApplyNonpositivePolicy(fit$nNonPositive, fit$n, nonpositive)

  fitted_ncp <- if (identical(target, "ncp")) fit$estimate else 0
  fitted_df <- if (identical(target, "ncp")) comparison$df else fit$estimate
  seed <- bootSeed %||% .ppeDefaultSeed(x, comparison)

  cvm <- .ppeCramerVonMises(fit$retained, df = fitted_df, ncp = fitted_ncp, target = target)
  p_value <- .ppeDiagnosticPValue(
    fit$retained, df = comparison$df, estimate = fit$estimate, target = target,
    bootstrapSamples = bootstrapSamples, seed = seed
  )

  retained_sorted <- sort(fit$retained)
  n <- fit$nRetained
  grid <- seq(0, max(retained_sorted) * 1.05, length.out = gridLength)
  envelope <- .ppeEnvelopeBounds(
    n = n, df = fitted_df, ncp = fitted_ncp, grid = grid,
    envelopeSamples = bootstrapSamples, conf.level = conf.level, seed = seed
  )

  # `df` here is the comparison's DECLARED/nominal value (the chi-square
  # reference used for the critical value), kept comparable across power and
  # type1 rows -- see ppeSummary()'s identical convention. It is NOT the
  # diagnosed value for a Type-I comparison: there, df is what was
  # ESTIMATED, and `df` above (comparison$df) can differ from it (e.g. a
  # nominal df = 1 with a fitted df = 1.0074). `parameter`/`estimate` (also
  # ppeSummary()'s own column names) are what the CvM statistic and p-value
  # are actually computed against, and must be reported alongside `df` so
  # the diagnosed value is never invisible in its own output.
  row <- data.frame(
    comparison = comparison$label, parameter = fit$parameter, estimate = fit$estimate,
    cvm = cvm, p_value = p_value,
    n = fit$n, n_nonpositive = fit$nNonPositive,
    df = comparison$df %||% NA_real_, df_source = comparison$dfSource,
    stringsAsFactors = FALSE
  )
  ecdf_data <- data.frame(
    x = retained_sorted, comparison = comparison$label, stringsAsFactors = FALSE
  )
  envelope_data <- data.frame(
    x = grid, lower = envelope$lower, upper = envelope$upper,
    fitted = stats::pchisq(grid, df = fitted_df, ncp = fitted_ncp),
    comparison = comparison$label, stringsAsFactors = FALSE
  )
  probs <- (seq_len(n) - 0.5) / n
  qq_data <- data.frame(
    theoretical = stats::qchisq(probs, df = fitted_df, ncp = fitted_ncp),
    empirical = retained_sorted,
    comparison = comparison$label,
    stringsAsFactors = FALSE
  )

  list(row = row, ecdf_data = ecdf_data, envelope_data = envelope_data, qq_data = qq_data)
}

# One subtitle line summarising the diagnostic table. A single comparison
# gets its own n/excluded/df detail (this is what the "the diagnostic
# subtitle names the excluded count" test checks); with several comparisons
# a per-comparison breakdown would not fit one subtitle line, so it
# aggregates instead -- the full per-comparison detail is always available
# from the "ppeDiagnostics" attribute, not just the subtitle text.
.ppeDiagnosticsSubtitle <- function(diagTable) {
  if (nrow(diagTable) == 1L) {
    row <- diagTable[1L, ]
    sprintf(
      "n = %d, %d excluded (nonpositive), df = %s (%s)",
      row$n, row$n_nonpositive, format(row$df), row$df_source %||% "NA"
    )
  } else {
    sprintf(
      "%d comparisons; %d of %d test statistics excluded (nonpositive) in total",
      nrow(diagTable), sum(diagTable$n_nonpositive), sum(diagTable$n)
    )
  }
}

#' Diagnose whether the fitted noncentral chi-square describes the observed test statistics
#'
#' `ppeSummary()`/`plotSSEPpePower()` treat the retained likelihood-ratio test
#' statistics as draws from a noncentral chi-square (Ueckert, Karlsson &
#' Hooker 2016; the method PsN implements as `distribution_mle`) and fit it
#' by maximum likelihood. Nothing in that pipeline checks whether the
#' noncentral chi-square actually describes the data -- this function does:
#' it computes a Cramer-von Mises discrepancy between the empirical CDF of
#' the retained statistics and the fitted chi-square CDF, with a
#' parametric-bootstrap p-value (draw synthetic data from the FITTED model,
#' refit, recompute the discrepancy, and see how often a correctly-specified
#' world produces one at least this large -- see `.ppeDiagnosticPValue()`),
#' and renders an ECDF-vs-fitted-CDF panel with a pointwise bootstrap
#' envelope alongside a QQ panel.
#'
#' Both `power` and `type1` comparisons can be diagnosed: a power
#' comparison's noncentrality is the value diagnosed (df held at the
#' comparison's known value); a Type-I comparison's df is the value
#' diagnosed (ncp held at 0) -- exactly the `target` split
#' `.ppeParametricBootstrap()` and `.ppeFit()` already use.
#'
#' This is diagnostic evidence about how well the noncentral chi-square
#' describes the observed data, not a pass/fail certification. A small
#' `p_value` says the fitted distribution is a poor description of the real
#' data; it does not by itself invalidate `ppeSummary()`'s point estimate,
#' and a large `p_value` does not prove the assumption exactly right, only
#' not detectably wrong at this sample size. Read it the way `.ppeFit()`'s
#' boundary-solution warning is read: a flag to weigh alongside the estimate,
#' not a gate on using it.
#'
#' @param x An `nlmixr2SSE` object.
#' @param comparisons Optional `sseComparison()` object or list of them. When
#'   `NULL`, uses `x$runInfo$comparisons`, then the legacy
#'   simulation-vs-alternative comparisons (with a warning that `df` was
#'   inferred rather than declared; see `.resolveComparisons()`). With more
#'   than one comparison, the plot facets by comparison and the
#'   `"ppeDiagnostics"` attribute has one row per comparison.
#' @param models Optional subset of alternative-model labels, used only when
#'   `comparisons` is `NULL`. Mutually exclusive with `comparisons`.
#' @param nonpositive What to do when a comparison has non-positive test
#'   statistics excluded from the fit: `"warn"` (the default), `"error"`, or
#'   `"drop"` (silent). See `.ppeApplyNonpositivePolicy()`.
#' @param conf.level Confidence level for the pointwise ECDF envelope.
#' @param bootstrapSamples Number of parametric-bootstrap replicates per
#'   comparison, used both for the p-value's null distribution and for the
#'   ECDF envelope.
#' @param bootSeed Optional integer seed for the bootstrap. `NULL` (the
#'   default) derives a seed deterministically from the run and each
#'   comparison's label (`.ppeDefaultSeed()`), so repeated calls on the same
#'   run reproduce the same diagnostic without the caller managing seeds.
#' @param ... Reserved for future extensions.
#'
#' @return A `ggplot` object (a `patchwork` combination of an ECDF panel and
#'   a QQ panel), with the numeric diagnostics attached as a
#'   `"ppeDiagnostics"` attribute (one row per comparison, with `comparison`,
#'   `parameter` (`"ncp"` or `"df"`, whichever was diagnosed), `estimate`
#'   (the fitted value of that parameter -- the value the discrepancy and
#'   p-value are actually about), `cvm`, `p_value`, `n`, `n_nonpositive`,
#'   `df` (the comparison's declared/nominal df, kept comparable across
#'   power and type1 rows; for a Type-I comparison this can differ from
#'   `estimate`), and `df_source`) so results are auditable without reading
#'   pixels.
#' @export
plotSSEPpeDiagnostics <- function(
  x,
  comparisons = NULL,
  models = NULL,
  nonpositive = c("warn", "error", "drop"),
  conf.level = 0.95,
  bootstrapSamples = 1000L,
  bootSeed = NULL,
  ...
) {
  .assertNamespace("ggplot2", "construct SSE plots")
  .assertNamespace("patchwork", "construct the PPE diagnostics panel")
  .assertSSEObject(x)
  nonpositive <- match.arg(nonpositive)

  # Resolved EXACTLY ONCE, with ppe = TRUE, so the "df inferred" warning
  # fires exactly once when it applies -- see .ppePowerPlotDataMle()'s header
  # in R/plot-sse.R for why re-resolving an already-resolved list would
  # silently swallow the warning rather than firing it twice.
  cmps <- .resolveComparisons(x, comparisons %||% x$runInfo$comparisons, models, ppe = TRUE)

  ineligible <- Filter(function(cmp) !isTRUE(cmp$ppeEligible), cmps)
  if (length(ineligible) > 0L) {
    bad <- vapply(ineligible, `[[`, character(1), "label")
    # cli::cli_abort() directly, not .abortSSE(): the "i" elements below are
    # already named, and .abortSSE() unconditionally wraps its argument in
    # c("!" = ...), turning an already-named "i" element into "!.i" -- a name
    # cli does not recognise as a bullet type. See .ppeFit()'s identical
    # workaround in R/sse-ppe.R.
    cli::cli_abort(c(
      "!" = "{length(bad)} comparison{?s} not eligible for distribution-based PPE: {.val {bad}}.",
      "i" = "PPE assumes a noncentral chi-square alternative, which an explicit {.arg criticalValue} gives no basis for.",
      "i" = "Build the comparison with {.arg df} instead of {.arg criticalValue}."
    ))
  }

  panels <- lapply(cmps, function(cmp) {
    .ppeDiagnosticsForComparison(
      x, cmp, nonpositive = nonpositive, conf.level = conf.level,
      bootstrapSamples = bootstrapSamples, bootSeed = bootSeed
    )
  })

  diag_table <- do.call(rbind, lapply(panels, `[[`, "row"))
  ecdf_data <- do.call(rbind, lapply(panels, `[[`, "ecdf_data"))
  envelope_data <- do.call(rbind, lapply(panels, `[[`, "envelope_data"))
  qq_data <- do.call(rbind, lapply(panels, `[[`, "qq_data"))

  facet_multi <- length(cmps) > 1L

  p_ecdf <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(
      data = envelope_data,
      ggplot2::aes(x = x, ymin = lower, ymax = upper, group = comparison),
      alpha = 0.2, fill = "grey40"
    ) +
    ggplot2::geom_line(
      data = envelope_data,
      ggplot2::aes(x = x, y = fitted, colour = comparison)
    ) +
    ggplot2::stat_ecdf(
      data = ecdf_data,
      ggplot2::aes(x = x, colour = comparison),
      geom = "step"
    ) +
    ggplot2::labs(
      x = "Test statistic", y = "Cumulative probability",
      colour = "Comparison", title = "ECDF vs fitted CDF"
    ) +
    ggplot2::theme_minimal()
  if (facet_multi) {
    p_ecdf <- p_ecdf + ggplot2::facet_wrap(~comparison)
  }

  p_qq <- ggplot2::ggplot(
    qq_data,
    ggplot2::aes(x = theoretical, y = empirical, colour = comparison)
  ) +
    ggplot2::geom_point(size = 1.2, alpha = 0.7) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey40") +
    ggplot2::labs(
      x = "Theoretical quantile", y = "Empirical quantile",
      colour = "Comparison", title = "QQ plot"
    ) +
    ggplot2::theme_minimal()
  if (facet_multi) {
    p_qq <- p_qq + ggplot2::facet_wrap(~comparison, scales = "free")
  }

  combined <- (p_ecdf | p_qq) + ggplot2::labs(
    subtitle = .ppeDiagnosticsSubtitle(diag_table),
    caption = paste(
      "The bootstrap envelope is pointwise, not simultaneous: it covers the",
      "fitted model's sampling variability at each point separately, not the",
      "whole curve at once."
    )
  )

  attr(combined, "ppeDiagnostics") <- diag_table
  combined
}
