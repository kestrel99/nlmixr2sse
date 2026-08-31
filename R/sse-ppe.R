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
#' `P(X > 0) = 1` under the noncentral chi-square, so this is a selection,
#' not a truncation: there is no missing normalization constant to restore.
#' Every discarded value is one the fitted family says cannot occur, so
#' retaining only positives makes this a selected-subset fit; its bias, if
#' any, depends on why those values occurred and is not knowably signed by
#' this fact alone. This matches the published method (PsN); the mitigation
#' is disclosure, not correction -- `nNonPositive` is always reported
#' alongside `estimate`, see `.ppeApplyNonpositivePolicy()`.
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
  # to the retained values. P(X > 0) = 1 under that model, so f(x | X > 0) =
  # f(x) exactly -- there is no truncation-renormalization issue. The
  # discarded values are ones the model says cannot occur, so retaining only
  # positives is a selection, not a truncation correction; the discarded
  # count must always be reported alongside the fit (see
  # .ppeApplyNonpositivePolicy()) so the exclusion is auditable.
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
    "i" = "This makes the fit a selected-subset estimator; its bias is not knowably signed by this fact alone."
  )
  switch(policy,
    error = cli::cli_abort(msg),
    warn = cli::cli_warn(msg),
    drop = invisible(NULL)
  )
}

# Parametric-bootstrap uncertainty for the distribution_mle point estimate.
#
# The MLE (.ppeChiSquareMle()) gives a point estimate only. This section adds
# an interval by resampling FROM THE FITTED MODEL: draw synthetic test
# statistics from the fitted noncentral chi-square, refit the MLE to each
# draw, and take percentile quantiles of the refit estimates.
#
# This is explicitly `interval_type = "model_based"`, never "empirical":
# `rchisq()` draws are strictly positive, so no replicate can ever reproduce
# the truncation the REAL data underwent (retaining only test statistics
# `> 0`; see `.ppeChiSquareMle()`'s header). The interval therefore covers
# estimator variability under the fitted model only -- never model
# misspecification, and never the truncation bias `.ppeChiSquareMle()`
# already documents as mild and upward. Keep this distinction in naming,
# comments, and user-facing docs; do not let "model_based" drift into
# "empirical" anywhere downstream.

#' Set the RNG seed for one expression, then restore the caller's stream exactly
#'
#' Mirrors the save/restore pattern already used by the `ppe_dofv()` test
#' helper: on entry, records whether `.Random.seed` exists and its value; on
#' exit, restores that exact state (including its absence, if it was absent)
#' regardless of how `expr` completes. `seed = NULL` skips seeding entirely
#' and runs `expr` against whatever stream the caller already has.
#'
#' @param seed A single integer, or `NULL` to leave the RNG untouched.
#' @param expr An expression to evaluate under the seeded stream.
#' @return The value of `expr`.
#' @noRd
.withPpeSeed <- function(seed, expr) {
  if (is.null(seed)) {
    return(force(expr))
  }
  has_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old <- if (has_seed) get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit({
    if (has_seed) {
      assign(".Random.seed", old, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  }, add = TRUE)
  set.seed(seed)
  force(expr)
}

#' Parametric-bootstrap interval for a distribution_mle point estimate
#'
#' Draws `bootstrapSamples` synthetic test-statistic samples from the FITTED
#' noncentral chi-square (`rchisq(nRetained, df, ncp = estimate)` when
#' `target = "ncp"`; `rchisq(nRetained, df = estimate, ncp = 0)` when
#' `target = "df"`), refits `.ppeChiSquareMle()` to each draw, and returns
#' percentile quantiles of the refit estimates. A refit that errors (e.g. a
#' convergence failure) is counted in `n_failed` rather than aborting the
#' whole bootstrap.
#'
#' @param estimate The fitted point estimate from `.ppeChiSquareMle()`
#'   (`ncp` when `target = "ncp"`, `df` when `target = "df"`).
#' @param df Degrees of freedom, used only when `target = "ncp"`.
#' @param nRetained Number of retained (positive) test statistics the
#'   original fit used; each bootstrap draw is the same size.
#' @param bootstrapSamples Number of bootstrap replicates. `0` skips the
#'   bootstrap entirely and returns `NA` bounds with `n_successful = 0`.
#' @param target Which parameter `estimate` is (`"ncp"` or `"df"`).
#' @param conf.level Confidence level for the percentile interval.
#' @param seed Passed to `.withPpeSeed()`; `NULL` leaves the caller's RNG
#'   stream untouched (still consumed, just not reset first) -- callers that
#'   need reproducibility supply an explicit seed.
#' @return A list with `ci_lower`, `ci_upper`, `n_successful`, `n_failed`,
#'   `interval_type` (always `"model_based"`, see the module note above), and
#'   `draws` (the successful refit estimates).
#' @noRd
.ppeParametricBootstrap <- function(estimate, df, nRetained,
                                    bootstrapSamples = 1000L,
                                    target = c("ncp", "df"),
                                    conf.level = 0.95, seed = NULL) {
  target <- match.arg(target)
  empty <- list(
    ci_lower = NA_real_, ci_upper = NA_real_, n_successful = 0L,
    n_failed = 0L, interval_type = "model_based", draws = numeric(0)
  )
  if (bootstrapSamples <= 0L) {
    return(empty)
  }

  draws <- .withPpeSeed(seed, {
    vapply(seq_len(bootstrapSamples), function(i) {
      tryCatch({
        if (identical(target, "ncp")) {
          .ppeChiSquareMle(
            stats::rchisq(nRetained, df = df, ncp = estimate),
            df = df
          )$estimate
        } else {
          .ppeChiSquareMle(
            stats::rchisq(nRetained, df = estimate, ncp = 0),
            ncp = 0
          )$estimate
        }
      }, error = function(e) NA_real_)
    }, numeric(1))
  })

  ok <- draws[is.finite(draws)]
  if (length(ok) == 0L) {
    return(empty)
  }
  a <- 1 - conf.level
  q <- stats::quantile(ok, probs = c(a / 2, 1 - a / 2), names = FALSE)
  list(
    ci_lower = q[[1L]],
    ci_upper = q[[2L]],
    n_successful = length(ok),
    n_failed = length(draws) - length(ok),
    interval_type = "model_based",
    draws = ok
  )
}

#' A per-comparison bootstrap seed, derived deterministically from the run
#'
#' Derived from the run's own seed (never truly random) so a given run plots
#' identically every time it is re-rendered, while remaining distinct across
#' runs and across comparisons within one run (via the label offset).
#'
#' The label offset is a POSITION-WEIGHTED sum of character codes, not a
#' plain `sum(utf8ToInt(label))`: a plain sum is permutation-invariant, so
#' two labels that are anagrams of each other collide on the same offset --
#' and not just in contrived cases. `"dose A vs B"` and `"dose B vs A"`, a
#' natural pair of labels for the two directions of one comparison, both sum
#' to the same total; weighting each character's code by its 1-based
#' position breaks that symmetry (verified: 887/887 under the plain sum
#' becomes 4870/4865 weighted). Colliding offsets mean colliding bootstrap
#' seeds -- and, when the two comparisons' point estimates also happen to
#' match, byte-identical bootstrap draws -- silently undermining the
#' "remaining distinct across comparisons" guarantee above.
#'
#' The addition is done in double precision and cast back to integer only at
#' the end: `runInfo$seed` is caller-supplied and can legitimately be close
#' to `.Machine$integer.max` (2147483647L). Doing `base + offset` in integer
#' arithmetic risks a silent overflow to `NA` (with a compiler warning R
#' surfaces as a real warning) rather than wrapping the way this function
#' intends; computing in `double` and reducing modulo 2147483647 before the
#' final `as.integer()` avoids that entirely, since the sum is always tiny
#' relative to double's ~15-digit precision.
#'
#' @param x An `nlmixr2SSE` object.
#' @param comparison A resolved `sseComparison` (must have `$label`).
#' @return A single integer seed.
#' @noRd
.ppeDefaultSeed <- function(x, comparison) {
  base <- as.double(x$runInfo$seed %||% 1L)
  codes <- utf8ToInt(comparison$label)
  offset <- sum(codes * seq_along(codes)) %% 100000L
  as.integer((base + offset) %% 2147483647)
}

#' Fit the distribution_mle estimator and its bootstrap interval for one comparison
#'
#' @param x An `nlmixr2SSE` object.
#' @param comparison A resolved `sseComparison`.
#' @param conf.level Confidence level for the bootstrap interval.
#' @param nonpositive Nonpositive-statistic policy, see
#'   `.ppeApplyNonpositivePolicy()`.
#' @param bootstrapSamples,bootSeed Passed to `.ppeParametricBootstrap()`.
#'   `bootSeed = NULL` derives a seed from the run via `.ppeDefaultSeed()`.
#' @return A list combining `.ppeChiSquareMle()`'s and
#'   `.ppeParametricBootstrap()`'s fields, plus `comparison`, `target`
#'   (`"ncp"` for a power comparison, `"df"` for a Type-I comparison), and
#'   `probability` (the estimated power/Type-I rate at the comparison's own
#'   `criticalValue`).
#' @noRd
.ppeFit <- function(x, comparison, conf.level = 0.95, nonpositive = "warn",
                    bootstrapSamples = 1000L, bootSeed = NULL) {
  if (!isTRUE(comparison$ppeEligible)) {
    # A pre-named vector (the "i" element below) passed through .abortSSE()
    # would double-wrap: .abortSSE() unconditionally does c("!" = ...), which
    # turns an already-named "i" element into "!.i" -- a name cli does not
    # recognise as a bullet type, silently losing that bullet's formatting.
    # Call cli::cli_abort() directly instead, as .ppeApplyNonpositivePolicy()
    # already does for the same reason.
    cli::cli_abort(c(
      "!" = "Comparison {.val {comparison$label}} has an explicit {.arg criticalValue} and no {.arg df}.",
      "i" = "Distribution-based PPE assumes a noncentral chi-square alternative, which a custom critical value gives no basis for."
    ))
  }
  stat <- .comparisonTestStatistic(x, comparison)
  target <- if (identical(comparison$mode, "type1")) "df" else "ncp"
  fit <- if (identical(target, "ncp")) {
    .ppeChiSquareMle(stat, df = comparison$df)
  } else {
    .ppeChiSquareMle(stat, ncp = 0)
  }
  .ppeApplyNonpositivePolicy(fit$nNonPositive, fit$n, nonpositive)
  if (isTRUE(fit$boundary)) {
    cli::cli_warn(c(
      "!" = "The {target} estimate for {.val {comparison$label}} is at its lower bound.",
      "i" = "Estimated power equals alpha and the interval is degenerate.",
      "i" = "This is the constrained maximum-likelihood estimate, not a numerical failure."
    ))
  }
  boot <- .ppeParametricBootstrap(
    fit$estimate, df = comparison$df, nRetained = fit$nRetained,
    bootstrapSamples = bootstrapSamples, target = target,
    conf.level = conf.level, seed = bootSeed %||% .ppeDefaultSeed(x, comparison)
  )
  probability <- if (identical(target, "ncp")) {
    stats::pchisq(comparison$criticalValue, df = comparison$df,
                  ncp = fit$estimate, lower.tail = FALSE)
  } else {
    stats::pchisq(comparison$criticalValue, df = fit$estimate,
                  ncp = 0, lower.tail = FALSE)
  }
  c(fit, boot, list(comparison = comparison, target = target,
                    probability = probability))
}

#' Model-based parametric power estimation summary
#'
#' One row per comparison, combining the `distribution_mle` point estimate
#' (`.ppeChiSquareMle()`) with its parametric-bootstrap interval
#' (`.ppeParametricBootstrap()`, always `interval_type = "model_based"`: the
#' bootstrap draws from the FITTED noncentral chi-square, so the interval
#' covers estimator variability under the fitted model only, never model
#' misspecification, and is never "empirical") and provenance (`df_source`,
#' `boundary`, retained/discarded counts). This is the summary-table
#' counterpart to `plotSSEPpePower()`'s curve/point-range rendering: same
#' estimator, same bootstrap, tabular rather than plotted.
#'
#' @param x An `nlmixr2SSE` object.
#' @param comparisons Optional `sseComparison()` object or list of them. When
#'   `NULL`, uses `x$runInfo$comparisons`, then the legacy
#'   simulation-vs-alternative comparisons (with a warning that `df` was
#'   inferred rather than declared; see `.resolveComparisons()`).
#' @param models Optional subset of alternative-model labels, used only when
#'   `comparisons` is `NULL`. Mutually exclusive with `comparisons`.
#' @param conf.level Confidence level for the bootstrap interval.
#' @param nonpositive What to do when a comparison has non-positive test
#'   statistics excluded from the fit: `"warn"` (the default), `"error"`, or
#'   `"drop"` (silent). See `.ppeApplyNonpositivePolicy()`.
#' @param bootstrapSamples Number of parametric-bootstrap replicates per
#'   comparison. `0` skips the bootstrap and returns `NA` intervals while
#'   still reporting the point estimate and all counts -- the cheap option
#'   when the interval is not needed.
#' @param bootSeed Optional integer seed for the bootstrap. `NULL` (the
#'   default) derives a seed deterministically from the run and the
#'   comparison's label (`.ppeDefaultSeed()`), so repeated calls on the same
#'   run reproduce the same interval without the caller managing seeds.
#' @return A `data.frame` with one row per comparison.
#' @export
ppeSummary <- function(x, comparisons = NULL, models = NULL,
                       conf.level = 0.95,
                       nonpositive = c("warn", "error", "drop"),
                       bootstrapSamples = 1000L, bootSeed = NULL) {
  .assertSSEObject(x)
  nonpositive <- match.arg(nonpositive)
  cmps <- .resolveComparisons(x, comparisons %||% x$runInfo$comparisons,
                              models, ppe = TRUE)
  rows <- lapply(cmps, function(cmp) {
    f <- .ppeFit(x, cmp, conf.level = conf.level, nonpositive = nonpositive,
                 bootstrapSamples = bootstrapSamples, bootSeed = bootSeed)
    data.frame(
      comparison = cmp$label, full = cmp$full, reduced = cmp$reduced,
      mode = cmp$mode, parameter = f$parameter,
      df = cmp$df %||% NA_real_, df_source = cmp$dfSource,
      alpha = cmp$alpha, critical_value = cmp$criticalValue,
      n = f$n, n_nonpositive = f$nNonPositive,
      estimate = f$estimate, ci_lower = f$ci_lower, ci_upper = f$ci_upper,
      interval_type = f$interval_type,
      n_bootstrap_successful = f$n_successful,
      n_bootstrap_failed = f$n_failed,
      boundary = f$boundary, probability = f$probability,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}
