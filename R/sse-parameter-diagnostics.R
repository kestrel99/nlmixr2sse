# Adequacy diagnostics for `parameterSource = "covariance"` replicate draws.
#
# `.resolveCovarianceParameterSets()` (R/sse-helpers.R) draws replicate-level
# THETA/OMEGA parameter sets from a local Normal/inverse-Wishart approximation
# around a real model fit. Both draw modes have known, DOCUMENTED artefacts on
# the raw parameter scale:
#
#  * `covarianceDraw = "joint"` transforms OMEGA to a log-Cholesky scale before
#    drawing multivariate-Normal (see R/sse-omega-joint.R), which is unbiased
#    on the TRANSFORMED scale but introduces upward mean drift on the RAW
#    scale: for a lognormal-shaped diagonal element, `E[exp(X)] =
#    exp(mu + sigma^2/2)`, i.e. the raw-scale mean is inflated by a factor of
#    `exp(se^2 / (2 * omega^2))` relative to the fitted value.
#  * `covarianceDraw = "independent_iw"` draws OMEGA from an inverse-Wishart
#    with a single `nu` per block (see `.omegaWishartSpec()`,
#    R/sse-omega-draw.R), so only the "binding" element -- whichever implies
#    the smallest `nu` -- matches its reported SE exactly; every other element
#    in the block comes out MORE dispersed than reported.
#
# `parameterDrawSummary()` surfaces both artefacts numerically, so a user
# auditing a covariance-mode run can see them instead of silently trusting the
# approximation.

# A realized mean more than this fraction away from its target mean is
# flagged. 5% is comfortably below the ~11-13% drift the "joint" mode's
# log-Cholesky back-transformation produces at a typical SE/Omega ratio of
# 0.5 (see the module note above), while still being well outside ordinary
# Monte Carlo noise at realistic replicate counts (a few hundred to a few
# thousand).
.parameterDrawMeanDriftThreshold <- 0.05

#' Summarize how well drawn covariance-mode parameters match their targets
#'
#' For a completed `runSSE()` run with `parameterSource = "covariance"`, this
#' reports -- for every generating parameter, across all replicates -- how the
#' ACTUALLY DRAWN values (`x$initialValues`) compare to the INTENDED targets
#' that seeded the draw (the fitted value and reported standard error from
#' `x$runInfo$parameterSourceInfo$targets`). This is a diagnostic only: it
#' never modifies, truncates, or resamples the recorded draws.
#'
#' Two known artefacts, both intentional consequences of the draw machinery
#' rather than bugs, are what this table is mainly for surfacing:
#'
#' * Under `covarianceDraw = "joint"`, OMEGA is drawn on a log-Cholesky scale
#'   (unconstrained, so a plain multivariate-Normal draw always back-
#'   transforms to a positive-definite matrix). That is unbiased on the
#'   transformed scale but inflates the RAW-scale mean upward -- watch
#'   `realized_mean` versus `target_mean` and `mean_drift_flag`.
#' * Under `covarianceDraw = "independent_iw"`, each OMEGA block is drawn from
#'   an inverse-Wishart with a single degrees-of-freedom `nu` (`binding_nu`),
#'   chosen so that exactly one diagonal element -- the "binding" one --
#'   matches its reported SE exactly. Every other element in the same block
#'   comes out more dispersed than reported -- watch `dispersion_ratio`
#'   (`realized_sd / target_sd`), which will exceed 1 for non-binding
#'   elements.
#'
#' THETA out-of-bound draws (`n_out_of_domain`) and OMEGA positive-
#' definiteness failures (`n_not_positive_definite`) are also counted per
#' parameter; both draw modes are positive-definite by construction (see
#' R/sse-omega-draw.R and R/sse-omega-joint.R), so the latter should be `0` in
#' virtually every realistic case, and a nonzero count is itself diagnostic of
#' an extreme, ill-conditioned fit.
#'
#' @param x An `nlmixr2SSE` object produced by a `parameterSource =
#'   "covariance"` `runSSE()` call.
#' @param ... Unused; reserved for future extension.
#' @return A `data.frame` with one row per generating parameter (THETA or
#'   OMEGA entry), with columns:
#'   \describe{
#'     \item{parameter}{The parameter name, as it appears in
#'       `x$initialValues$parameter`.}
#'     \item{n}{Finite draw count.}
#'     \item{realized_mean, realized_sd, realized_median, realized_q025,
#'       realized_q975, realized_min, realized_max}{Summary statistics of the
#'       actually-drawn values across replicates.}
#'     \item{target_mean, target_sd}{The Normal/inverse-Wishart target mean
#'       and standard deviation the draw was seeded from, when known.}
#'     \item{dispersion_ratio}{`realized_sd / target_sd`; noticeably above 1
#'       for `independent_iw`'s non-binding OMEGA elements.}
#'     \item{mean_drift_flag}{`TRUE` when `realized_mean` deviates from
#'       `target_mean` by more than 5% relatively (a fixed internal
#'       threshold, comfortably below the joint mode's typical drift and
#'       comfortably above ordinary Monte Carlo noise).}
#'     \item{lower, upper}{Recoverable finite THETA bounds from the model's
#'       `ini()` block (`-Inf`/`Inf` when unconstrained); `NA` for OMEGA
#'       entries, which have no such domain.}
#'     \item{n_out_of_domain}{Count of replicates whose THETA draw fell
#'       outside `[lower, upper]`. Diagnostic only -- draws are never
#'       truncated or resampled. `NA` for parameters with no tracked bounds.}
#'     \item{binding_nu}{The `independent_iw` inverse-Wishart degrees of
#'       freedom for this parameter's OMEGA block (shared by every parameter
#'       in the same block). `NA` for THETA parameters and for OMEGA blocks
#'       with no usable reported SE.}
#'     \item{n_not_positive_definite}{Count of replicates where this
#'       parameter's OMEGA block failed a positive-definiteness check
#'       (shared by every parameter in the same block). `NA` for THETA
#'       parameters.}
#'   }
#'
#'   For `covarianceDraw = "joint"` runs with two or more drawn parameters,
#'   the result also carries a `"jointDependence"` attribute: a `data.frame`
#'   of pairwise empirical correlations between the drawn parameters, computed
#'   from `x$initialValues`. This is a descriptive fact about what actually
#'   came out of the draw, not a claim that it reproduces `fit$cov` exactly.
#' @export
parameterDrawSummary <- function(x, ...) {
  .assertSSEObject(x)
  if (!identical(x$runInfo$parameterSource, "covariance")) {
    .abortSSE(
      "{.arg x} must be a {.val covariance}-parameterSource run; got {.val {x$runInfo$parameterSource %||% 'NA'}}."
    )
  }
  .parameterDrawSummaryCore(x$runInfo, x$initialValues)
}

#' The actual work behind `parameterDrawSummary()`, taking `runInfo` and
#' `initialValues` directly rather than a full `nlmixr2SSE` object.
#'
#' Split out so `runSSE()` can precompute and cache this table for a
#' covariance-mode run (see `.newNlmixr2SSE()`'s `parameterDrawSummary`
#' field) without needing to assemble a complete `nlmixr2SSE` object first.
#' @noRd
.parameterDrawSummaryCore <- function(runInfo, initialValues) {
  long <- initialValues
  if (is.null(long) || nrow(long) == 0L) {
    .abortSSE("{.arg x$initialValues} has no recorded replicate draws to summarize.")
  }

  targets <- runInfo$parameterSourceInfo$targets
  parameters <- unique(long$parameter)

  rows <- lapply(parameters, function(p) {
    values <- long$value[long$parameter == p]
    values <- values[is.finite(values)]
    .parameterDrawSummaryRow(p, values, targets)
  })
  summary_tab <- do.call(rbind, rows)
  summary_tab <- summary_tab[order(summary_tab$parameter), , drop = FALSE]
  rownames(summary_tab) <- NULL

  # Positive-definiteness and binding-nu are BLOCK-level facts (every
  # parameter belonging to the same OMEGA block shares one value), inferred
  # from the "omega(<row>,<col>)" naming convention alone so this works even
  # when parameterSourceInfo carries no explicit block structure (as in a
  # minimal, hand-built fixture).
  summary_tab <- .annotateOmegaBlocks(summary_tab, long, targets)

  if (identical(runInfo$covarianceDraw, "joint")) {
    attr(summary_tab, "jointDependence") <- .parameterDrawJointDependence(long)
  }

  summary_tab
}

#' One parameter's realized-vs-target summary row
#' @noRd
.parameterDrawSummaryRow <- function(parameter, values, targets) {
  target_row <- if (!is.null(targets)) {
    targets[targets$parameter == parameter, , drop = FALSE]
  } else {
    NULL
  }
  get_target <- function(col, default = NA_real_) {
    if (is.null(target_row) || nrow(target_row) == 0L || !col %in% names(target_row)) {
      return(default)
    }
    target_row[[col]][[1L]]
  }

  target_mean <- get_target("target_mean")
  target_sd <- get_target("target_sd")
  lower <- get_target("lower", NA_real_)
  # `upper` is not part of the minimal target schema every caller is expected
  # to populate; its absence means "not tracked", which is the same as
  # "unbounded above" for out-of-domain counting -- not a gap to work around.
  upper <- get_target("upper", Inf)
  if (is.na(upper)) {
    upper <- Inf
  }

  realized_mean <- .safeMean(values)
  realized_sd <- .safeSd(values)

  mean_drift_flag <- if (is.na(target_mean) || is.na(realized_mean) ||
                            target_mean == 0) {
    NA
  } else {
    abs(realized_mean - target_mean) / abs(target_mean) >
      .parameterDrawMeanDriftThreshold
  }

  # Bounds are "tracked" whenever `lower`/`upper` carry an actual value --
  # even an unconstrained theta (lower = -Inf, upper = Inf) is legitimately
  # checked and legitimately never out-of-domain, so that counts as tracked.
  # NA on both sides means "no domain to check" (e.g. an OMEGA entry), and
  # yields NA_integer_ rather than a misleading 0.
  bounds_tracked <- !is.na(lower) || is.finite(upper)
  n_out_of_domain <- if (!bounds_tracked) {
    NA_integer_
  } else {
    lo <- if (is.na(lower)) -Inf else lower
    as.integer(sum(values < lo | values > upper))
  }

  data.frame(
    parameter = parameter,
    n = length(values),
    realized_mean = realized_mean,
    realized_sd = realized_sd,
    realized_median = .safeMedian(values),
    realized_q025 = if (length(values) == 0L) NA_real_ else stats::quantile(values, 0.025, names = FALSE, type = 7),
    realized_q975 = if (length(values) == 0L) NA_real_ else stats::quantile(values, 0.975, names = FALSE, type = 7),
    realized_min = .safeMin(values),
    realized_max = .safeMax(values),
    target_mean = target_mean,
    target_sd = target_sd,
    dispersion_ratio = if (is.na(target_sd) || target_sd == 0) NA_real_ else realized_sd / target_sd,
    mean_drift_flag = mean_drift_flag,
    lower = lower,
    upper = upper,
    n_out_of_domain = n_out_of_domain,
    binding_nu = get_target("binding_nu", NA_real_),
    n_not_positive_definite = NA_integer_,
    stringsAsFactors = FALSE
  )
}

#' Group `omega(<row>,<col>)` parameter names into independent blocks
#'
#' Union-find over eta names read straight off the parameter labels
#' themselves, mirroring `.omegaBlocks()` (R/sse-omega-draw.R) but operating
#' on names rather than `fit$ui` indices -- so it needs nothing beyond
#' `x$initialValues`, and works even for a hand-built fixture that carries no
#' block structure of its own.
#' @return named list, each element a character vector of covNames in that block
#' @noRd
.parameterDrawOmegaGroups <- function(omegaParams) {
  m <- regmatches(
    omegaParams,
    regexec("^omega\\(([^,]+),([^,]+)\\)$", omegaParams)
  )
  rowEta <- vapply(m, function(z) if (length(z) == 3L) z[[2L]] else NA_character_, character(1))
  colEta <- vapply(m, function(z) if (length(z) == 3L) z[[3L]] else NA_character_, character(1))
  keep <- !is.na(rowEta) & !is.na(colEta)

  etas <- unique(c(rowEta[keep], colEta[keep]))
  parent <- stats::setNames(etas, etas)
  findRoot <- function(e) {
    while (parent[[e]] != e) {
      e <- parent[[e]]
    }
    e
  }
  for (i in which(keep)) {
    a <- findRoot(rowEta[[i]])
    b <- findRoot(colEta[[i]])
    if (a != b) {
      parent[[a]] <- b
    }
  }

  roots <- vapply(etas, findRoot, character(1))
  # Every off-diagonal entry already unions its row and column eta above, so
  # a parameter's row-eta root alone identifies its block.
  unname(split(omegaParams[keep], roots[rowEta[keep]]))
}

#' Attach per-block `binding_nu`/`n_not_positive_definite` to the summary table
#' @noRd
.annotateOmegaBlocks <- function(summaryTab, long, targets) {
  omega_params <- grep("^omega\\(", summaryTab$parameter, value = TRUE)
  if (length(omega_params) == 0L) {
    return(summaryTab)
  }
  groups <- .parameterDrawOmegaGroups(omega_params)

  for (grp in groups) {
    pd_count <- .omegaBlockNotPositiveDefiniteCount(grp, long)
    summaryTab$n_not_positive_definite[summaryTab$parameter %in% grp] <- pd_count

    # binding_nu: if targets already carries it, per-parameter values were
    # honoured in .parameterDrawSummaryRow(); this only fills gaps by
    # propagating whichever non-NA value the block already has, so a target
    # table that annotates only one representative row still yields a
    # complete column.
    have <- summaryTab$binding_nu[summaryTab$parameter %in% grp]
    non_na <- have[!is.na(have)]
    if (length(non_na) > 0L) {
      summaryTab$binding_nu[summaryTab$parameter %in% grp] <- non_na[[1L]]
    }
  }

  summaryTab
}

#' Count replicates where an OMEGA block's drawn matrix is not positive-definite
#'
#' Reconstructs each replicate's block matrix from the long-format draws and
#' checks `eigen(..., only.values = TRUE)$values > 0`. Both draw modes are
#' positive-definite by construction (see the module note in
#' R/sse-omega-draw.R and R/sse-omega-joint.R), so this is expected to be `0`
#' in virtually every realistic case -- but it is checked for real rather
#' than assumed, since extreme condition numbers can violate that in floating
#' point.
#' @noRd
.omegaBlockNotPositiveDefiniteCount <- function(covNames, long) {
  m <- regmatches(covNames, regexec("^omega\\(([^,]+),([^,]+)\\)$", covNames))
  rowEta <- vapply(m, `[[`, character(1), 2L)
  colEta <- vapply(m, `[[`, character(1), 3L)
  etas <- unique(c(rowEta, colEta))
  p <- length(etas)

  sub <- long[long$parameter %in% covNames, , drop = FALSE]
  wide <- stats::reshape(
    sub[, c("replicate", "parameter", "value")],
    idvar = "replicate", timevar = "parameter", direction = "wide"
  )
  names(wide) <- sub("^value\\.", "", names(wide))

  if (p == 1L) {
    v <- wide[[covNames[[1L]]]]
    return(sum(!is.na(v) & v <= 0))
  }

  bad <- 0L
  for (i in seq_len(nrow(wide))) {
    mat <- matrix(NA_real_, p, p, dimnames = list(etas, etas))
    ok <- TRUE
    for (k in seq_along(covNames)) {
      val <- wide[[covNames[[k]]]][[i]]
      if (is.na(val)) {
        ok <- FALSE
        break
      }
      mat[rowEta[[k]], colEta[[k]]] <- val
      mat[colEta[[k]], rowEta[[k]]] <- val
    }
    if (!ok || anyNA(mat)) {
      next
    }
    ev <- suppressWarnings(eigen(mat, symmetric = TRUE, only.values = TRUE)$values)
    if (!all(is.finite(ev)) || min(ev) <= 0) {
      bad <- bad + 1L
    }
  }
  bad
}

#' Empirical pairwise correlations between jointly-drawn parameters
#'
#' Descriptive only: a "joint" draw incorporates `fit$cov`'s THETA<->OMEGA
#' covariance through a nonlinear (log-Cholesky) transform, so the realized
#' RAW-scale correlation is not expected to reproduce `fit$cov` exactly, even
#' though the transformed-scale draw is exact multivariate-Normal.
#' @noRd
.parameterDrawJointDependence <- function(long) {
  params <- unique(long$parameter)
  if (length(params) < 2L) {
    return(data.frame(
      parameter1 = character(0), parameter2 = character(0),
      correlation = numeric(0), stringsAsFactors = FALSE
    ))
  }
  wide <- stats::reshape(
    long[, c("replicate", "parameter", "value")],
    idvar = "replicate", timevar = "parameter", direction = "wide"
  )
  names(wide) <- sub("^value\\.", "", names(wide))

  pairs <- utils::combn(params, 2L, simplify = FALSE)
  rows <- lapply(pairs, function(pr) {
    x <- wide[[pr[[1L]]]]
    y <- wide[[pr[[2L]]]]
    ok <- is.finite(x) & is.finite(y)
    corr <- if (sum(ok) < 2L) NA_real_ else stats::cor(x[ok], y[ok])
    data.frame(
      parameter1 = pr[[1L]], parameter2 = pr[[2L]],
      correlation = corr, stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}
