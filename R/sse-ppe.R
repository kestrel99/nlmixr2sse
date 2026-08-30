# Maximum-likelihood parametric power estimation (PPE).
#
# Ueckert, Karlsson & Hooker (2016) -- the method PsN implements -- treat the
# likelihood-ratio test statistics from an SSE run as draws from a
# noncentral chi-square with known degrees of freedom, and fit ONE
# noncentrality parameter by maximum likelihood from the whole distribution,
# then scale it linearly with study size. This is the statistical core of
# that estimator; see `method = "distribution_mle"` in `.ppePowerPlotData()`
# (R/plot-sse.R) for how it is wired into the power-vs-sample-size curve.
#
# `method = "exceedance"` remains available: for each threshold separately it
# counts how many statistics exceed it and inverts that single exceedance
# probability to get a noncentrality. That uses one bit of information per
# threshold and can imply a different effect size at every threshold -- the
# defining weakness the MLE fixes.

#' Fit a noncentral chi-square to retained test statistics by maximum likelihood
#'
#' Retains only test statistics `> 0` (the noncentral chi-square has no
#' support at or below zero) and fits either the noncentrality parameter
#' (`ncp`, with `df` known) or `df` (with `ncp` known, used to sanity-check
#' recovery) by direct maximisation of the unconditional noncentral
#' chi-square log-likelihood.
#'
#' The fit deliberately does NOT renormalise for the truncation: since
#' `P(X > 0) = 1` under the noncentral chi-square, every discarded value is
#' one the model says cannot occur, so the estimate is mildly biased upward.
#' This matches the published method (PsN); the mitigation is disclosure, not
#' correction -- `nNonPositive` is always reported alongside `estimate`, see
#' `.ppeApplyNonpositivePolicy()`.
#'
#' A retained mean below `df` is a legitimate outcome: the constrained MLE
#' sits at its lower bound (`boundary = TRUE`, `estimate` effectively `0`,
#' i.e. estimated power equals alpha). This is the constrained MLE behaving
#' correctly, not a numerical failure -- it is flagged, never corrected.
#'
#' @param dofv Numeric vector of test statistics (`OFV_reduced - OFV_full`).
#'   Non-finite values are dropped before counting; of the finite values,
#'   those `<= 0` are excluded from the fit and counted in `nNonPositive`.
#' @param df Degrees of freedom, when `ncp` is to be estimated. Exactly one
#'   of `df`/`ncp` must be supplied; the other is estimated.
#' @param ncp Noncentrality, when `df` is to be estimated instead.
#' @return A list with `estimate`, `parameter` (`"ncp"` or `"df"`, whichever
#'   was estimated), `objective` (negative log-likelihood at the optimum),
#'   `convergence`, `message`, `retained` (the values the fit used), `n`
#'   (finite input values), `nRetained`, `nNonPositive`, and `boundary`
#'   (`TRUE` when the optimum sits at the lower bound).
#' @noRd
.ppeChiSquareMle <- function(dofv, df = NULL, ncp = NULL) {
  estimate_ncp <- is.null(ncp)
  # Validate the df/ncp contract BEFORE counting retained values. The plan's
  # original ordering checked "at least 2 positive values" first, so
  # `.ppeChiSquareMle(c(-1, 1))` (neither df nor ncp supplied, and only 1
  # positive value) reported "needs at least 2 positive test statistics" --
  # true but not the actual problem. A caller who forgot both arguments
  # should be told that.
  if (estimate_ncp == is.null(df)) {
    .abortSSE("Supply exactly one of {.arg df} and {.arg ncp}; the other is estimated.")
  }

  finite <- dofv[is.finite(dofv)]
  retained <- finite[finite > 0]
  n_nonpositive <- length(finite) - length(retained)
  if (length(retained) < 2L) {
    .abortSSE(c(
      "Distribution-based PPE needs at least 2 positive test statistics.",
      "i" = "Got {length(retained)} positive of {length(finite)} finite values.",
      "i" = "If all values are negative, check that the comparison names the models the right way round."
    ))
  }

  # PsN and Ueckert (2016) fit the unconditional noncentral chi-square density
  # to the retained values. P(X > 0) = 1 under that model, so the discarded
  # values are ones it says cannot occur; the fit is therefore mildly biased
  # upward and the discarded count must always be reported alongside it.
  objective <- if (estimate_ncp) {
    function(par) -sum(stats::dchisq(retained, df = df, ncp = par, log = TRUE))
  } else {
    function(par) -sum(stats::dchisq(retained, df = par, ncp = ncp, log = TRUE))
  }
  init <- mean(retained) - if (estimate_ncp) df else ncp
  fit <- stats::optim(
    par = init, fn = objective, lower = 1e-16,
    method = "L-BFGS-B"
  )
  if (!identical(as.integer(fit$convergence), 0L)) {
    .abortSSE(c(
      "The test-statistic likelihood did not converge.",
      "i" = "{.field optim} said: {fit$message %||% 'no message'}."
    ))
  }
  list(
    estimate = fit$par,
    parameter = if (estimate_ncp) "ncp" else "df",
    objective = fit$value,
    convergence = as.integer(fit$convergence),
    message = fit$message %||% NA_character_,
    retained = retained,
    n = length(finite),
    nRetained = length(retained),
    nNonPositive = n_nonpositive,
    boundary = fit$par <= 1e-8
  )
}

#' Warn about, abort on, or silently accept non-positive test statistics
#'
#' Never touches the counts themselves -- `.ppeChiSquareMle()` always
#' computes and reports `nNonPositive`/`nRetained`; this only controls what
#' happens when that count is nonzero.
#'
#' @param nNonPositive,n Counts from `.ppeChiSquareMle()` (`nNonPositive` and
#'   `n`, the number of finite input values).
#' @param policy One of `"warn"`, `"error"`, `"drop"`.
#' @return `invisible(NULL)`.
#' @noRd
.ppeApplyNonpositivePolicy <- function(nNonPositive, n, policy) {
  if (nNonPositive == 0L) {
    return(invisible(NULL))
  }

  # cli's `{?s}`/`{?is/are}` pluralisation keys off the LAST numeric
  # interpolation before the marker by default. With both `nNonPositive` and
  # `n` interpolated before it, that reads off `n` (the denominator), not
  # `nNonPositive` -- e.g. nNonPositive = 1, n = 5 would render "1 of 5 test
  # statistics are not positive", plural agreeing with the wrong number.
  # cli::qty() pins the pluralisation to `nNonPositive` explicitly; verified
  # by hand that (1, 5) then reads "1 of 5 test statistic is not positive".
  #
  # Built with the leading "!" name already present (unlike .abortSSE()'s
  # single-string call sites, which rely on it to inject that name): passing
  # this vector THROUGH .abortSSE() would double-wrap it -- .abortSSE() does
  # `c("!" = ...)` unconditionally, turning an already-named "!" element into
  # "!.!" and "i" into "!.i", names cli does not recognise as bullet types,
  # so the rendered message silently loses its bullets. Calling cli::cli_abort()
  # directly here avoids that; both branches use exactly the same message.
  msg <- c(
    "!" = paste0(
      "{nNonPositive} of {n} {cli::qty(nNonPositive)}test statistic{?s} ",
      "{cli::qty(nNonPositive)}{?is/are} not positive."
    ),
    "i" = "The noncentral chi-square has no support there, so they are excluded from the fit.",
    "i" = "The estimate is biased upward when this fraction is large."
  )
  switch(policy,
    error = cli::cli_abort(msg),
    warn = cli::cli_warn(msg),
    drop = invisible(NULL)
  )
}
