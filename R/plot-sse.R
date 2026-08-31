#' Plot SSE results
#'
#' `plot()` dispatches to one of the package's Phase 7 plotting helpers.
#'
#' @param x An `nlmixr2SSE` object.
#' @param type Plot type. Supported values are `"parameter_bias"`,
#'   `"parameter_estimates"`, `"ofv_distribution"`, `"power"`,
#'   `"ppe_power"`, and `"diagnostics"`.
#' @param ... Additional arguments passed through to the selected plot helper.
#'
#' @return A `ggplot` object, or a patchwork object for diagnostics.
#' @name plot-nlmixr2SSE
NULL

.assertSSEObject <- function(x, arg = "x") {
  if (!inherits(x, "nlmixr2SSE")) {
    .abortSSE(paste0("`", arg, "` must inherit from `nlmixr2SSE`."))
  }
  invisible(x)
}

.assertNamespace <- function(pkg, reason) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    .abortSSE(
      paste0("Package `", pkg, "` is required to ", reason, ".")
    )
  }
  invisible(TRUE)
}

.filterPlotModels <- function(data, models = NULL) {
  if (is.null(models)) {
    return(data)
  }
  checkmate::assertCharacter(models, any.missing = FALSE, min.len = 1L)
  data[data$model_label %in% models, , drop = FALSE]
}

.filterPlotParameters <- function(data, parameters = NULL) {
  if (is.null(parameters)) {
    return(data)
  }
  checkmate::assertCharacter(parameters, any.missing = FALSE, min.len = 1L)
  data[data$parameter %in% parameters, , drop = FALSE]
}

.parameterBiasPlotData <- function(
  x,
  statistic = "relative_bias",
  models = NULL,
  parameters = NULL,
  includeUnmatched = FALSE
) {
  .assertSSEObject(x)
  checkmate::assertChoice(
    statistic,
    c("bias", "relative_bias", "relative_absolute_bias", "rmse", "relative_rmse")
  )
  checkmate::assertFlag(includeUnmatched)

  data <- x$parameterSummary %||% .emptyParameterSummary()
  data <- data[data$statistic == statistic, , drop = FALSE]
  data <- .filterPlotModels(data, models = models)
  data <- .filterPlotParameters(data, parameters = parameters)
  if (!isTRUE(includeUnmatched)) {
    data <- data[data$matched, , drop = FALSE]
  }
  data[order(data$parameter, data$model_label), , drop = FALSE]
}

.parameterEstimatePlotData <- function(
  x,
  models = NULL,
  parameters = NULL,
  includeUnmatched = FALSE
) {
  .assertSSEObject(x)
  checkmate::assertFlag(includeUnmatched)

  header <- attr(x$rawResults, "rawResultsHeader", exact = TRUE)
  if (is.null(header)) {
    .abortSSE("Raw results are missing canonical header metadata.")
  }

  param_cols <- c(
    header$theta_cols %||% character(0),
    header$omega_cols %||% character(0),
    header$sigma_cols %||% character(0)
  )
  if (length(param_cols) == 0L) {
    return(data.frame(
      sample = integer(0),
      model_label = character(0),
      role = character(0),
      parameter = character(0),
      estimate = numeric(0),
      truth = numeric(0),
      matched = logical(0),
      stringsAsFactors = FALSE
    ))
  }

  row_mask <- .summaryMask(
    x$rawResults,
    outFilter = x$runInfo$control$outFilter %||% NULL
  )
  raw <- x$rawResults[x$rawResults$sample > 0L & row_mask, , drop = FALSE]

  spec_map <- .specMapFromSnapshots(x$runInfo$fitSpecs %||% list())
  truths <- x$initialValues
  rows <- vector("list", nrow(raw) * length(param_cols))
  idx <- 1L

  for (i in seq_len(nrow(raw))) {
    model_label <- raw$model_label[[i]]
    model_params <- names(.schemaParamPresence(spec_map[[model_label]] %||% list(
      thetaCols = character(0),
      omegaCols = character(0),
      sigmaCols = character(0)
    )))
    for (param in param_cols) {
      matched <- param %in% model_params
      if (!isTRUE(includeUnmatched) && !matched) {
        next
      }
      truth_row <- truths$replicate == raw$sample[[i]] & truths$parameter == param
      rows[[idx]] <- data.frame(
        sample = as.integer(raw$sample[[i]]),
        model_label = model_label,
        role = raw$role[[i]],
        parameter = param,
        estimate = as.numeric(raw[[param]][[i]]),
        truth = if (any(truth_row)) truths$value[match(TRUE, truth_row)] else NA_real_,
        matched = matched,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }

  data <- if (idx == 1L) {
    rows[0]
  } else {
    do.call(rbind, rows[seq_len(idx - 1L)])
  }
  data <- .filterPlotModels(data, models = models)
  data <- .filterPlotParameters(data, parameters = parameters)
  data[order(data$parameter, data$model_label, data$sample), , drop = FALSE]
}

.ofvDeltaPlotData <- function(x, models = NULL) {
  .assertSSEObject(x)

  fit_specs <- x$runInfo$fitSpecs %||% list()
  alternative_specs <- Filter(
    function(spec) identical(spec$role, "alternative"),
    fit_specs
  )
  if (length(alternative_specs) == 0L) {
    return(data.frame(
      sample = integer(0),
      model_label = character(0),
      delta_ofv = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  samples <- sort(unique(x$rawResults$sample[x$rawResults$sample > 0L]))
  row_mask <- .summaryMask(
    x$rawResults,
    outFilter = x$runInfo$control$outFilter %||% NULL
  )
  filtered <- x$rawResults[row_mask, , drop = FALSE]

  ofv_by_label <- lapply(fit_specs, function(spec) {
    rows <- filtered$sample > 0L & filtered$model_label == spec$label
    vec <- rep(NA_real_, length(samples))
    if (any(rows)) {
      vec[match(filtered$sample[rows], samples)] <- filtered$objf[rows]
    }
    vec
  })
  names(ofv_by_label) <- vapply(fit_specs, `[[`, character(1), "label")

  ref_ofv <- if (!is.null(x$runInfo$control$refOfv)) {
    rep(as.numeric(x$runInfo$control$refOfv), length(samples))
  } else if (isTRUE(x$runInfo$estimateSimulation)) {
    sim_label <- fit_specs[[match(
      "simulation",
      vapply(fit_specs, `[[`, character(1), "role")
    )]]$label
    ofv_by_label[[sim_label]]
  } else {
    ofv_by_label[[alternative_specs[[1L]]$label]]
  }

  rows <- vector("list", length(alternative_specs))
  for (i in seq_along(alternative_specs)) {
    label <- alternative_specs[[i]]$label
    rows[[i]] <- data.frame(
      sample = samples,
      model_label = label,
      delta_ofv = ref_ofv - ofv_by_label[[label]],
      stringsAsFactors = FALSE
    )
  }

  data <- do.call(rbind, rows)
  data <- data[is.finite(data$delta_ofv), , drop = FALSE]
  data <- .filterPlotModels(data, models = models)
  data[order(data$model_label, data$sample), , drop = FALSE]
}

.powerPlotData <- function(x, thresholds = NULL, direction = "power", models = NULL) {
  .assertSSEObject(x)
  checkmate::assertChoice(direction, c("power", "type1"))

  data <- x$ofvSummary
  data <- data[!is.na(data$direction) & data$direction == direction, , drop = FALSE]
  if (!is.null(thresholds)) {
    checkmate::assertNumeric(thresholds, any.missing = FALSE, min.len = 1L)
    data <- data[data$threshold %in% thresholds, , drop = FALSE]
  }
  data <- .filterPlotModels(data, models = models)
  data[order(data$model_label, data$threshold), , drop = FALSE]
}

.studySampleSize <- function(x) {
  .assertSSEObject(x)

  size <- x$runInfo$studySampleSize %||% NA_integer_
  if (!is.finite(size) || size < 1L) {
    .abortSSE(
      "PPE plotting requires {.field runInfo$studySampleSize} metadata from a completed {.fn runSSE} call."
    )
  }

  list(
    size = as.integer(size),
    unit = x$runInfo$studySampleUnit %||% "subjects"
  )
}

.simulationSpec <- function(x) {
  fit_specs <- x$runInfo$fitSpecs %||% list()
  sim_idx <- match(
    "simulation",
    vapply(fit_specs, `[[`, character(1), "role")
  )
  if (is.na(sim_idx)) {
    .abortSSE("Could not locate the simulation-model specification in `runInfo$fitSpecs`.")
  }
  fit_specs[[sim_idx]]
}

.specParameterCount <- function(spec) {
  schema <- spec$schema %||% list()
  length(unique(c(
    schema$thetaCols %||% character(0),
    schema$omegaCols %||% character(0),
    schema$sigmaCols %||% character(0)
  )))
}

.modelDegreesFreedom <- function(x, model_label) {
  fit_specs <- x$runInfo$fitSpecs %||% list()
  spec_idx <- match(model_label, vapply(fit_specs, `[[`, character(1), "label"))
  if (is.na(spec_idx)) {
    return(1L)
  }

  df <- .specParameterCount(.simulationSpec(x)) - .specParameterCount(fit_specs[[spec_idx]])
  as.integer(max(df, 1L))
}

.ppeTailProbability <- function(threshold, df, ncp) {
  stats::pchisq(threshold, df = df, ncp = ncp, lower.tail = FALSE)
}

.clipProbability <- function(prob, n) {
  eps <- 0.5 / (n + 1)
  min(max(prob, eps), 1 - eps)
}

.binomialInterval <- function(successes, total, conf.level = 0.95) {
  alpha <- 1 - conf.level
  lower <- if (successes <= 0L) {
    0
  } else {
    stats::qbeta(alpha / 2, successes, total - successes + 1)
  }
  upper <- if (successes >= total) {
    1
  } else {
    stats::qbeta(1 - alpha / 2, successes + 1, total - successes)
  }
  c(lower = lower, upper = upper)
}

.solveNcpForProbability <- function(probability, threshold, df) {
  base_prob <- .ppeTailProbability(threshold, df = df, ncp = 0)
  if (probability <= base_prob) {
    return(0)
  }

  upper <- 1
  upper_prob <- .ppeTailProbability(threshold, df = df, ncp = upper)
  while (upper_prob < probability && upper < 1e6) {
    upper <- upper * 2
    upper_prob <- .ppeTailProbability(threshold, df = df, ncp = upper)
  }
  if (upper_prob < probability) {
    return(upper)
  }

  stats::uniroot(
    function(ncp) {
      .ppeTailProbability(threshold, df = df, ncp = ncp) - probability
    },
    interval = c(0, upper),
    tol = 1e-8
  )$root
}

.defaultPpeStudySizes <- function(baseStudySize, pointNcp, threshold, df, targetPower) {
  target_prob <- targetPower / 100
  end_size <- baseStudySize

  if (is.finite(pointNcp) && pointNcp > 0 && target_prob > 0 && target_prob < 1) {
    target_ncp <- .solveNcpForProbability(
      probability = target_prob,
      threshold = threshold,
      df = df
    )
    if (is.finite(target_ncp) && target_ncp > 0) {
      end_size <- max(
        baseStudySize,
        ceiling(baseStudySize * target_ncp / pointNcp)
      )
    }
  }

  if (end_size <= 300L) {
    seq_len(end_size)
  } else {
    unique(as.integer(round(seq(1, end_size, length.out = 200L))))
  }
}

.emptyPpePowerPlotData <- function(method) {
  base <- data.frame(
    model_label = character(0),
    threshold = numeric(0),
    degrees_freedom = integer(0),
    base_study_size = integer(0),
    study_size = integer(0),
    power = numeric(0),
    power_lower = numeric(0),
    power_upper = numeric(0),
    stringsAsFactors = FALSE
  )
  if (identical(method, "exceedance")) {
    base$threshold_exceedance_probability <- numeric(0)
    base$threshold_implied_ncp <- numeric(0)
  } else {
    base$mle_ncp <- numeric(0)
    base$mle_n_retained <- integer(0)
    base$mle_n_nonpositive <- integer(0)
    base$mle_boundary <- logical(0)
    base$df_source <- character(0)
  }
  base
}

# Today's estimator, unchanged: for each threshold separately, invert that
# threshold's own exceedance rate into a noncentrality. Kept byte-identical
# under `method = "exceedance"` -- the only addition is exposing the
# per-threshold probability/noncentrality that were already computed
# internally, under names that can never be confused with the single
# whole-distribution noncentrality `method = "distribution_mle"` fits.
.ppePowerPlotDataExceedance <- function(
  x, delta_data, thresholds, base_study, targetPower, conf.level, studySizes
) {
  rows <- list()
  row_index <- 1L
  for (model_label in unique(delta_data$model_label)) {
    model_rows <- delta_data$model_label == model_label
    test_stat <- -delta_data$delta_ofv[model_rows]
    test_stat <- test_stat[is.finite(test_stat)]
    if (length(test_stat) == 0L) {
      next
    }

    df <- .modelDegreesFreedom(x, model_label)
    total <- length(test_stat)

    for (threshold in thresholds) {
      successes <- sum(test_stat > threshold)
      point_prob <- .clipProbability(successes / total, total)
      point_ncp <- .solveNcpForProbability(point_prob, threshold = threshold, df = df)
      interval <- .binomialInterval(successes, total, conf.level = conf.level)
      lower_ncp <- .solveNcpForProbability(interval[["lower"]], threshold = threshold, df = df)
      upper_ncp <- .solveNcpForProbability(interval[["upper"]], threshold = threshold, df = df)

      sizes <- studySizes %||% .defaultPpeStudySizes(
        baseStudySize = base_study$size,
        pointNcp = point_ncp,
        threshold = threshold,
        df = df,
        targetPower = targetPower
      )

      scaled_point <- point_ncp * sizes / base_study$size
      scaled_lower <- lower_ncp * sizes / base_study$size
      scaled_upper <- upper_ncp * sizes / base_study$size

      rows[[row_index]] <- data.frame(
        model_label = model_label,
        threshold = threshold,
        degrees_freedom = df,
        base_study_size = base_study$size,
        study_size = sizes,
        power = 100 * vapply(
          scaled_point,
          .ppeTailProbability,
          numeric(1),
          threshold = threshold,
          df = df
        ),
        power_lower = 100 * vapply(
          scaled_lower,
          .ppeTailProbability,
          numeric(1),
          threshold = threshold,
          df = df
        ),
        power_upper = 100 * vapply(
          scaled_upper,
          .ppeTailProbability,
          numeric(1),
          threshold = threshold,
          df = df
        ),
        threshold_exceedance_probability = point_prob,
        threshold_implied_ncp = point_ncp,
        stringsAsFactors = FALSE
      )
      row_index <- row_index + 1L
    }
  }

  if (row_index == 1L) {
    return(.emptyPpePowerPlotData("exceedance"))
  }

  data <- do.call(rbind, rows)
  data[order(data$threshold, data$model_label, data$study_size), , drop = FALSE]
}

# The new default: one noncentrality per comparison, fit by maximum
# likelihood from the WHOLE retained distribution (.ppeChiSquareMle()), then
# scaled linearly with study size. Unlike the exceedance estimator, this fit
# does not depend on `threshold` at all -- it is computed once per comparison
# and reused for every threshold, which is the entire point of the change.
#
# Restricted to power-mode comparisons: a Type-I comparison's estimated
# quantity is df, not a noncentrality that scales with study size, so it has
# no sample-size curve to draw here at all. plotSSEPpePower() renders any
# Type-I comparisons separately as a point-range
# (.ppeType1PointRangeData()/.plotSSEPpeType1PointRange()); a Type-I row that
# reaches this function (e.g. because a caller invoked it directly) is
# silently skipped rather than forced into a curve that means nothing for it.
.ppePowerPlotDataMle <- function(
  x, comparisons, models, thresholds, base_study, targetPower, studySizes,
  nonpositivePolicy, conf.level = 0.95, bootstrapSamples = 1000L, bootSeed = NULL,
  resolvedComparisons = NULL
) {
  # ppe = TRUE: distribution_mle treats df as KNOWN -- it is the exponent in
  # the noncentral chi-square density the whole fit rests on, not a cosmetic
  # label. Silencing this (as an earlier version of this function did, to
  # match the exceedance branch's long-standing silence) would leave no
  # signal anywhere -- no warning, no output column -- that df was guessed
  # from parameter counts rather than declared, for the method that is now
  # the package default. Exceedance's silence is the legacy posture this
  # whole feature replaces, not a precedent to match. `dfSource` is also
  # copied into the `df_source` output column below so provenance survives
  # even when the caller suppresses the warning.
  #
  # resolvedComparisons, when supplied, is used AS-IS instead of resolving
  # `comparisons %||% x$runInfo$comparisons` again. plotSSEPpePower() passes
  # its own already-resolved, mode-split list here: re-resolving a non-NULL
  # list would never re-enter .resolveComparisons()'s comparisons-is-NULL
  # branch, which is the ONLY place .legacyComparisons() (and its
  # inferred-df warning) ever runs -- so a second resolve of an
  # already-resolved list would silently make the warning never fire at all,
  # not fire it a second time. Bypassing resolution here (rather than
  # resolving twice and suppressing the second) is what keeps the warning
  # firing exactly once, from wherever the FIRST resolve happens.
  cmps <- resolvedComparisons %||% .resolveComparisons(
    x, comparisons %||% x$runInfo$comparisons, models = models, ppe = TRUE
  )

  ineligible <- Filter(function(cmp) !isTRUE(cmp$ppeEligible), cmps)
  if (length(ineligible) > 0L) {
    bad <- vapply(ineligible, `[[`, character(1), "label")
    .abortSSE(c(
      "{length(bad)} comparison{?s} not eligible for distribution-based PPE: {.val {bad}}.",
      "i" = "PPE assumes a noncentral chi-square alternative, which an explicit {.arg criticalValue} gives no basis for.",
      "i" = "Build the comparison with {.arg df} instead of {.arg criticalValue}, or use {.arg method = \"exceedance\"}."
    ))
  }

  # Type-I comparisons have no sample-size curve (df is not a noncentrality
  # that scales with N) -- see the function header. Skip them here; any
  # eligibility problem they have was already reported above regardless of
  # mode, so this filter only ever drops rows that were always going to be
  # rendered elsewhere, never ones that silently vanish unexplained.
  cmps <- Filter(function(cmp) identical(cmp$mode, "power"), cmps)

  rows <- list()
  row_index <- 1L
  for (cmp in cmps) {
    test_stat <- .comparisonTestStatistic(x, cmp)
    if (length(test_stat) == 0L) {
      next
    }

    mle_fit <- .ppeChiSquareMle(test_stat, df = cmp$df)
    .ppeApplyNonpositivePolicy(mle_fit$nNonPositive, mle_fit$n, nonpositivePolicy)

    # One bootstrap per comparison, reused across every threshold below --
    # the fit itself does not depend on threshold, so refitting per-threshold
    # would just repeat the same (expensive) work `length(thresholds)` times.
    boot <- .ppeParametricBootstrap(
      mle_fit$estimate, df = cmp$df, nRetained = mle_fit$nRetained,
      bootstrapSamples = bootstrapSamples, target = "ncp",
      conf.level = conf.level, seed = bootSeed %||% .ppeDefaultSeed(x, cmp)
    )
    has_interval <- is.finite(boot$ci_lower) && is.finite(boot$ci_upper)

    # mode is always "power" here (Type-I rows were filtered out above), so
    # the non-simulation member is always cmp$reduced -- mirrors what
    # .ofvDeltaPlotData()'s `model_label` always was for the legacy
    # power-mode comparisons the exceedance branch still builds.
    other_label <- cmp$reduced

    for (threshold in thresholds) {
      sizes <- studySizes %||% .defaultPpeStudySizes(
        baseStudySize = base_study$size,
        pointNcp = mle_fit$estimate,
        threshold = threshold,
        df = cmp$df,
        targetPower = targetPower
      )

      scale <- sizes / base_study$size
      scaled_ncp <- mle_fit$estimate * scale

      # The noncentrality scales linearly with study size, and pchisq(...,
      # ncp = ) is monotone increasing in ncp, so the bootstrap ncp bounds
      # map to the power bounds in the same order -- no reordering needed,
      # unlike the Type-I df bootstrap (.ppeType1PointRangeData()), where the
      # mapping is monotone DECREASING.
      power_lower <- rep(NA_real_, length(sizes))
      power_upper <- rep(NA_real_, length(sizes))
      if (has_interval) {
        power_lower <- 100 * vapply(
          boot$ci_lower * scale, .ppeTailProbability, numeric(1),
          threshold = threshold, df = cmp$df
        )
        power_upper <- 100 * vapply(
          boot$ci_upper * scale, .ppeTailProbability, numeric(1),
          threshold = threshold, df = cmp$df
        )
      }

      rows[[row_index]] <- data.frame(
        model_label = other_label,
        threshold = threshold,
        degrees_freedom = cmp$df,
        base_study_size = base_study$size,
        study_size = sizes,
        power = 100 * vapply(
          scaled_ncp,
          .ppeTailProbability,
          numeric(1),
          threshold = threshold,
          df = cmp$df
        ),
        power_lower = power_lower,
        power_upper = power_upper,
        mle_ncp = mle_fit$estimate,
        mle_n_retained = mle_fit$nRetained,
        mle_n_nonpositive = mle_fit$nNonPositive,
        mle_boundary = mle_fit$boundary,
        # "explicit" for a user-supplied sseComparison(df =), "parameter_count"
        # when df was inferred by .legacyComparisons() -- the same provenance
        # the ppe = TRUE warning above reports, kept in the data too so it
        # survives suppressWarnings()/expect_warning() at the call site.
        df_source = cmp$dfSource %||% NA_character_,
        stringsAsFactors = FALSE
      )
      row_index <- row_index + 1L
    }
  }

  if (row_index == 1L) {
    return(.emptyPpePowerPlotData("distribution_mle"))
  }

  data <- do.call(rbind, rows)
  data[order(data$threshold, data$model_label, data$study_size), , drop = FALSE]
}

# Type-I companion to .ppePowerPlotDataMle(): one row per Type-I comparison
# per threshold (NOT per study size -- df is not a noncentrality that scales
# with N, so there is no curve to grid over; see that function's header).
# The point estimate at threshold tau generalises
# ppeSummary()'s single-critical-value `probability` the same way the power
# curve generalises beyond a single threshold: P(T > tau) under the fitted
# df, evaluated at every tau in `thresholds`, exactly like power-mode's
# .ppeTailProbability() evaluates its scaled ncp at every tau rather than
# only at the comparison's own criticalValue.
.emptyPpeType1PlotData <- function() {
  data.frame(
    comparison = character(0),
    model_label = character(0),
    threshold = numeric(0),
    degrees_freedom = numeric(0),
    mle_df = numeric(0),
    rate = numeric(0),
    rate_lower = numeric(0),
    rate_upper = numeric(0),
    alpha = numeric(0),
    mle_n_retained = integer(0),
    mle_n_nonpositive = integer(0),
    mle_boundary = logical(0),
    df_source = character(0),
    stringsAsFactors = FALSE
  )
}

.ppeType1PointRangeData <- function(
  x, comparisons, thresholds, nonpositivePolicy, conf.level = 0.95,
  bootstrapSamples, bootSeed
) {
  rows <- list()
  row_index <- 1L
  for (cmp in comparisons) {
    test_stat <- .comparisonTestStatistic(x, cmp)
    if (length(test_stat) == 0L) {
      next
    }

    mle_fit <- .ppeChiSquareMle(test_stat, ncp = 0)
    .ppeApplyNonpositivePolicy(mle_fit$nNonPositive, mle_fit$n, nonpositivePolicy)

    boot <- .ppeParametricBootstrap(
      mle_fit$estimate, df = cmp$df, nRetained = mle_fit$nRetained,
      bootstrapSamples = bootstrapSamples, target = "df",
      conf.level = conf.level, seed = bootSeed %||% .ppeDefaultSeed(x, cmp)
    )
    has_interval <- is.finite(boot$ci_lower) && is.finite(boot$ci_upper)

    for (threshold in thresholds) {
      rate <- 100 * .ppeTailProbability(threshold, df = mle_fit$estimate, ncp = 0)

      rate_lower <- NA_real_
      rate_upper <- NA_real_
      if (has_interval) {
        # A larger df pushes the chi-square rightward, RAISING the tail
        # probability at a fixed threshold (measured: 0.019 at df = 0.5, up
        # to 1.000 at df = 32, holding threshold fixed) -- so the low end of
        # the df bootstrap interval maps to the low end of the rate
        # interval. min/max below does not rely on that direction, though:
        # it keeps rate_lower <= rate_upper regardless of which way the
        # mapping runs, defensive coding rather than a load-bearing claim
        # about monotonicity.
        r_lo_df <- 100 * .ppeTailProbability(threshold, df = boot$ci_lower, ncp = 0)
        r_hi_df <- 100 * .ppeTailProbability(threshold, df = boot$ci_upper, ncp = 0)
        rate_lower <- min(r_lo_df, r_hi_df)
        rate_upper <- max(r_lo_df, r_hi_df)
      }

      rows[[row_index]] <- data.frame(
        comparison = cmp$label,
        model_label = cmp$full,
        threshold = threshold,
        degrees_freedom = cmp$df,
        mle_df = mle_fit$estimate,
        rate = rate,
        rate_lower = rate_lower,
        rate_upper = rate_upper,
        alpha = cmp$alpha,
        mle_n_retained = mle_fit$nRetained,
        mle_n_nonpositive = mle_fit$nNonPositive,
        mle_boundary = mle_fit$boundary,
        df_source = cmp$dfSource %||% NA_character_,
        stringsAsFactors = FALSE
      )
      row_index <- row_index + 1L
    }
  }

  if (row_index == 1L) {
    return(.emptyPpeType1PlotData())
  }

  data <- do.call(rbind, rows)
  data[order(data$threshold, data$comparison), , drop = FALSE]
}

.ppePowerPlotData <- function(
  x,
  thresholds = NULL,
  models = NULL,
  studySizes = NULL,
  targetPower = 99,
  conf.level = 0.95,
  method = c("distribution_mle", "exceedance"),
  comparisons = NULL,
  nonpositivePolicy = c("warn", "error", "drop"),
  bootstrapSamples = 1000L,
  bootSeed = NULL,
  # Set only by plotSSEPpePower(), which resolves comparisons itself (once,
  # unsuppressed) to split them by mode before this function ever runs; see
  # .ppePowerPlotDataMle()'s header for why passing the resolved list through
  # here -- rather than back through `comparisons` -- is what keeps the
  # inferred-df warning firing exactly once instead of zero times.
  resolvedComparisons = NULL
) {
  .assertSSEObject(x)
  method <- match.arg(method)
  nonpositivePolicy <- match.arg(nonpositivePolicy)
  checkmate::assertNumber(targetPower, lower = 1, upper = 99.9, finite = TRUE)
  checkmate::assertNumber(conf.level, lower = 0.5, upper = 0.999, finite = TRUE)
  if (!is.null(studySizes)) {
    checkmate::assertIntegerish(
      studySizes,
      lower = 1,
      any.missing = FALSE,
      min.len = 1L
    )
    studySizes <- sort(unique(as.integer(studySizes)))
  }

  base_study <- .studySampleSize(x)
  delta_data <- .ofvDeltaPlotData(x, models = models)
  if (nrow(delta_data) == 0L) {
    return(.emptyPpePowerPlotData(method))
  }

  thresholds <- thresholds %||% sort(unique(x$powerSummary$threshold[x$powerSummary$threshold > 0]))
  checkmate::assertNumeric(thresholds, any.missing = FALSE, min.len = 1L, lower = 0)
  thresholds <- sort(unique(as.numeric(thresholds[thresholds > 0])))
  if (length(thresholds) == 0L) {
    .abortSSE("PPE plotting requires at least one positive OFV threshold.")
  }

  if (identical(method, "exceedance")) {
    return(.ppePowerPlotDataExceedance(
      x = x, delta_data = delta_data, thresholds = thresholds,
      base_study = base_study, targetPower = targetPower,
      conf.level = conf.level, studySizes = studySizes
    ))
  }

  .ppePowerPlotDataMle(
    x = x, comparisons = comparisons, models = models, thresholds = thresholds,
    base_study = base_study, targetPower = targetPower, studySizes = studySizes,
    nonpositivePolicy = nonpositivePolicy, conf.level = conf.level,
    bootstrapSamples = bootstrapSamples, bootSeed = bootSeed,
    resolvedComparisons = resolvedComparisons
  )
}

#' @rdname plot-nlmixr2SSE
#' @export
plot.nlmixr2SSE <- function(x, type = c(
  "parameter_bias",
  "parameter_estimates",
  "ofv_distribution",
  "power",
  "ppe_power",
  "diagnostics"
), ...) {
  type <- match.arg(type)

  switch(type,
    parameter_bias = plotSSEParameterBias(x, ...),
    parameter_estimates = plotSSEParameterEstimates(x, ...),
    ofv_distribution = plotSSEOfvDistribution(x, ...),
    power = plotSSEPower(x, ...),
    ppe_power = plotSSEPpePower(x, ...),
    diagnostics = plotSSEDiagnostics(x, ...)
  )
}

#' Plot SSE parameter bias summaries
#'
#' @param x An `nlmixr2SSE` object.
#' @param statistic Parameter-summary statistic to visualize.
#' @param models Optional model labels to retain.
#' @param parameters Optional parameter labels to retain.
#' @param includeUnmatched Logical. Include rows where a parameter is absent in
#'   a given hypothesis.
#' @param ... Reserved for future extensions.
#'
#' @return A `ggplot` object.
#' @export
plotSSEParameterBias <- function(
  x,
  statistic = "relative_bias",
  models = NULL,
  parameters = NULL,
  includeUnmatched = FALSE,
  ...
) {
  .assertNamespace("ggplot2", "construct SSE plots")
  .assertSSEObject(x)

  data <- .parameterBiasPlotData(
    x,
    statistic = statistic,
    models = models,
    parameters = parameters,
    includeUnmatched = includeUnmatched
  )

  if (nrow(data) == 0L) {
    return(
      ggplot2::ggplot() +
        ggplot2::labs(title = "SSE parameter summary", subtitle = "No data to display") +
        ggplot2::theme_void()
    )
  }

  ggplot2::ggplot(
    data,
    ggplot2::aes(x = model_label, y = value, fill = role)
  ) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::facet_wrap(~ parameter, scales = "free_y") +
    ggplot2::labs(
      x = "Model",
      y = statistic,
      fill = "Role",
      title = "SSE parameter summary"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
}

#' Plot SSE replicate-level parameter estimates
#'
#' @param x An `nlmixr2SSE` object.
#' @param models Optional model labels to retain.
#' @param parameters Optional parameter labels to retain.
#' @param includeTruth Logical. Overlay the simulated reference value for each
#'   replicate.
#' @param includeUnmatched Logical. Include rows where a parameter is absent in
#'   a given hypothesis.
#' @param ... Reserved for future extensions.
#'
#' @return A `ggplot` object.
#' @export
plotSSEParameterEstimates <- function(
  x,
  models = NULL,
  parameters = NULL,
  includeTruth = TRUE,
  includeUnmatched = FALSE,
  ...
) {
  .assertNamespace("ggplot2", "construct SSE plots")
  .assertSSEObject(x)
  checkmate::assertFlag(includeTruth)

  data <- .parameterEstimatePlotData(
    x,
    models = models,
    parameters = parameters,
    includeUnmatched = includeUnmatched
  )

  if (nrow(data) == 0L) {
    return(
      ggplot2::ggplot() +
        ggplot2::labs(title = "Replicate-level parameter estimates", subtitle = "No data to display") +
        ggplot2::theme_void()
    )
  }

  p <- ggplot2::ggplot(
    data,
    ggplot2::aes(x = model_label, y = estimate, fill = role)
  ) +
    ggplot2::geom_boxplot(outlier.alpha = 0.3) +
    ggplot2::geom_point(
      ggplot2::aes(group = model_label),
      position = ggplot2::position_jitter(width = 0.15, height = 0),
      size = 1.2,
      alpha = 0.5
    ) +
    ggplot2::facet_wrap(~ parameter, scales = "free_y") +
    ggplot2::labs(
      x = "Model",
      y = "Estimate",
      fill = "Role",
      title = "Replicate-level parameter estimates"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))

  if (isTRUE(includeTruth) && nrow(data) > 0L) {
    truth_data <- unique(data[, c("sample", "parameter", "truth"), drop = FALSE])
    p <- p + ggplot2::geom_point(
      data = truth_data,
      ggplot2::aes(x = 1, y = truth),
      inherit.aes = FALSE,
      shape = 4,
      size = 1.8,
      alpha = 0.7
    )
  }

  p
}

#' Plot SSE OFV distributions
#'
#' @param x An `nlmixr2SSE` object.
#' @param models Optional model labels to retain.
#' @param style Plot style: empirical CDF or histogram.
#' @param ... Reserved for future extensions.
#'
#' @return A `ggplot` object.
#' @export
plotSSEOfvDistribution <- function(
  x,
  models = NULL,
  style = c("ecdf", "histogram"),
  ...
) {
  .assertNamespace("ggplot2", "construct SSE plots")
  .assertSSEObject(x)
  style <- match.arg(style)

  data <- .ofvDeltaPlotData(x, models = models)

  p <- ggplot2::ggplot(
    data,
    ggplot2::aes(x = delta_ofv, colour = model_label, fill = model_label)
  ) +
    ggplot2::labs(
      x = "delta_OFV",
      y = if (identical(style, "ecdf")) "ECDF" else "Count",
      colour = "Model",
      fill = "Model",
      title = "SSE OFV distribution"
    ) +
    ggplot2::theme_minimal()

  if (identical(style, "ecdf")) {
    p + ggplot2::stat_ecdf(geom = "step", linewidth = 0.8)
  } else {
    p + ggplot2::geom_histogram(alpha = 0.35, position = "identity", bins = 20L)
  }
}

#' Plot SSE power or Type I curves
#'
#' @param x An `nlmixr2SSE` object.
#' @param thresholds Optional OFV thresholds to retain.
#' @param direction Which OFV direction to plot: `"power"` or `"type1"`.
#' @param models Optional model labels to retain.
#' @param ... Reserved for future extensions.
#'
#' @return A `ggplot` object.
#' @export
plotSSEPower <- function(
  x,
  thresholds = NULL,
  direction = c("power", "type1"),
  models = NULL,
  ...
) {
  .assertNamespace("ggplot2", "construct SSE plots")
  .assertSSEObject(x)
  direction <- match.arg(direction)

  data <- .powerPlotData(
    x,
    thresholds = thresholds,
    direction = direction,
    models = models
  )

  ggplot2::ggplot(
    data,
    ggplot2::aes(x = threshold, y = value, colour = model_label)
  ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::labs(
      x = "delta_OFV threshold",
      y = "Percent",
      colour = "Model",
      title = if (identical(direction, "power")) {
        "SSE power curve"
      } else {
        "SSE Type I curve"
      }
  ) +
    ggplot2::theme_minimal()
}

#' Plot a Type-I point-range panel for distribution_mle Type-I comparisons
#'
#' A Type-I comparison estimates `df` (with `ncp` held at `0`), not a
#' noncentrality -- there is nothing to scale with study size, so unlike
#' `.ppePowerPlotDataMle()`'s curve this renders one point-range per
#' comparison per threshold: the estimated exceedance probability under the
#' fitted `df` at that threshold (generalising `ppeSummary()`'s
#' single-critical-value `probability` across every threshold in
#' `thresholds`, exactly as the power curve generalises across thresholds
#' rather than using only the comparison's own `criticalValue`), with a
#' parametric-bootstrap point-range and a dashed reference line at each
#' comparison's own declared `alpha` (per-comparison, matching how
#' `comparisonSummary()`/`ppeSummary()` already treat `alpha`).
#'
#' @noRd
.plotSSEPpeType1PointRange <- function(
  x, comparisons, thresholds, nonpositivePolicy, conf.level = 0.95,
  bootstrapSamples, bootSeed
) {
  ineligible <- Filter(function(cmp) !isTRUE(cmp$ppeEligible), comparisons)
  if (length(ineligible) > 0L) {
    bad <- vapply(ineligible, `[[`, character(1), "label")
    # cli::cli_abort() directly, not .abortSSE(): the "i" elements below are
    # already named, and .abortSSE() unconditionally wraps its argument in
    # c("!" = ...), which turns an already-named "i" element into "!.i" --
    # a name cli does not recognise as a bullet type. See .ppeFit()'s same
    # workaround in R/sse-ppe.R.
    cli::cli_abort(c(
      "!" = "{length(bad)} comparison{?s} not eligible for distribution-based PPE: {.val {bad}}.",
      "i" = "PPE assumes a noncentral chi-square alternative, which an explicit {.arg criticalValue} gives no basis for.",
      "i" = "Build the comparison with {.arg df} instead of {.arg criticalValue}, or use {.arg method = \"exceedance\"}."
    ))
  }

  thresholds <- thresholds %||% sort(unique(x$powerSummary$threshold[x$powerSummary$threshold > 0]))
  checkmate::assertNumeric(thresholds, any.missing = FALSE, min.len = 1L, lower = 0)
  thresholds <- sort(unique(as.numeric(thresholds[thresholds > 0])))
  if (length(thresholds) == 0L) {
    .abortSSE("PPE plotting requires at least one positive OFV threshold.")
  }

  data <- .ppeType1PointRangeData(
    x, comparisons = comparisons, thresholds = thresholds,
    nonpositivePolicy = nonpositivePolicy, conf.level = conf.level,
    bootstrapSamples = bootstrapSamples, bootSeed = bootSeed
  )

  if (nrow(data) == 0L) {
    return(
      ggplot2::ggplot() +
        ggplot2::labs(title = "Parametric power estimation", subtitle = "No data to display") +
        ggplot2::theme_void()
    )
  }

  has_interval <- any(is.finite(data$rate_lower) & is.finite(data$rate_upper))

  p <- ggplot2::ggplot(
    data,
    ggplot2::aes(x = comparison, y = rate, colour = model_label)
  )
  p <- if (has_interval) {
    p + ggplot2::geom_pointrange(
      ggplot2::aes(ymin = rate_lower, ymax = rate_upper)
    )
  } else {
    p + ggplot2::geom_point(size = 2)
  }

  p +
    ggplot2::geom_hline(
      ggplot2::aes(yintercept = 100 * alpha),
      linetype = 2,
      linewidth = 0.5,
      colour = "grey40"
    ) +
    ggplot2::facet_wrap(
      ~ threshold,
      labeller = ggplot2::label_bquote(delta[OFV] == .(threshold))
    ) +
    ggplot2::labs(
      x = "Comparison",
      y = "Estimated Type-I rate (%)",
      colour = "Model",
      title = "Parametric power estimation (Type-I)"
    ) +
    ggplot2::coord_cartesian(ylim = c(0, 100)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
}

#' Plot PPE power-versus-sample-size curves
#'
#' @param x An `nlmixr2SSE` object.
#' @param thresholds Optional positive OFV thresholds to retain. By default,
#'   uses all positive thresholds from `powerSummary`.
#' @param models Optional model labels to retain.
#' @param studySizes Optional integer vector of study sample sizes to evaluate.
#'   When `NULL`, a default grid is generated from size `1` up to the estimated
#'   sample size achieving `targetPower`.
#' @param targetPower Target power percentage used when constructing the default
#'   `studySizes` grid.
#' @param conf.level Confidence level for the Monte-Carlo uncertainty ribbon
#'   (`"exceedance"`) or the parametric-bootstrap interval
#'   (`"distribution_mle"`).
#' @param method PPE estimator. `"distribution_mle"` (the default) fits ONE
#'   noncentrality parameter (power comparisons) or `df` (Type-I comparisons)
#'   to the whole retained test-statistic distribution by maximum likelihood
#'   (Ueckert, Karlsson & Hooker 2016; the method PsN implements). A power
#'   comparison's noncentrality is scaled linearly with study size into a
#'   curve; its `power_lower`/`power_upper` ribbon comes from a
#'   parametric-bootstrap interval (`interval_type = "model_based"`; see
#'   `.ppeParametricBootstrap()` -- it covers estimator variability under the
#'   fitted model only, never model misspecification, and is never
#'   "empirical"). A Type-I comparison has no sample-size curve -- `df` is not
#'   a noncentrality that scales -- so it renders separately as a
#'   point-range of the estimated Type-I rate against a dashed nominal
#'   `alpha` reference line; when a call mixes power and Type-I comparisons,
#'   the Type-I ones are omitted from the curve with a warning (see
#'   `ppeSummary()` for their estimates instead). `"exceedance"` is the
#'   original estimator: for each threshold separately, it inverts that
#'   threshold's own exceedance rate into a noncentrality, which can imply a
#'   different effect size at every threshold; it has no Type-I rendering and
#'   ignores comparison `mode` entirely.
#' @param comparisons Optional `sseComparison()` object or list of them,
#'   passed to [.resolveComparisons()][sseComparison()]. Only used under
#'   `method = "distribution_mle"`; defaults to `x$runInfo$comparisons`, then
#'   to the legacy simulation-vs-alternative comparisons.
#' @param nonpositivePolicy What to do when a comparison has non-positive test
#'   statistics excluded from the `"distribution_mle"` fit: `"warn"` (the
#'   default), `"error"`, or `"drop"` (silent). See
#'   `.ppeApplyNonpositivePolicy()`.
#' @param bootstrapSamples Number of parametric-bootstrap replicates per
#'   comparison, used only under `method = "distribution_mle"`. `0` skips the
#'   bootstrap and leaves the ribbon/point-range bounds `NA`.
#' @param bootSeed Optional integer seed for the bootstrap. `NULL` (the
#'   default) derives a seed deterministically from the run and each
#'   comparison's label (`.ppeDefaultSeed()`), so repeated calls on the same
#'   run reproduce the same interval without the caller managing seeds.
#' @param ... Reserved for future extensions.
#'
#' @section Limitations:
#' `method = "distribution_mle"` fits the UNCONDITIONAL noncentral
#' chi-square density to only the RETAINED (positive) test statistics,
#' excluding the rest -- it does not renormalize for that truncation, so the
#' fit is biased upward whenever non-positive statistics are common (the
#' excluded count is always reported so this is auditable, not hidden; see
#' `ppeSummary()`'s `n`/`n_nonpositive` columns).
#'
#' The `power_lower`/`power_upper` ribbon (and `ppeSummary()`'s bootstrap
#' interval) is a parametric bootstrap under the FITTED model
#' (`interval_type = "model_based"`): it covers estimator variability given
#' that the noncentral chi-square is correct, never the possibility that it
#' is not. It is not a statement about model misspecification.
#'
#' Extrapolating the curve to a study size that was never directly simulated
#' assumes the noncentrality scales linearly with size,
#' \eqn{\lambda(n) = \lambda_0 n / n_0}. This is an assumption, not something
#' this function verifies -- test it directly with `validateSSEPpeScaling()`
#' across runs simulated at more than one study size before trusting an
#' extrapolated value.
#'
#' Building an `sseComparison()` does not, by itself, prove the two named
#' models are nested or that their OFVs are on a comparable scale -- that is
#' the caller's responsibility. A hypothesis on a parameter boundary (e.g.
#' testing a variance against zero) does not follow ordinary chi-square
#' asymptotics; such comparisons need a user-supplied `criticalValue` (which
#' disables PPE, since PPE assumes a noncentral chi-square alternative) or an
#' empirical null calibrated by direct SSE, not the default `df`-based
#' reference.
#'
#' @return A `ggplot` object.
#' @export
plotSSEPpePower <- function(
  x,
  thresholds = NULL,
  models = NULL,
  studySizes = NULL,
  targetPower = 99,
  conf.level = 0.95,
  method = c("distribution_mle", "exceedance"),
  comparisons = NULL,
  nonpositivePolicy = c("warn", "error", "drop"),
  bootstrapSamples = 1000L,
  bootSeed = NULL,
  ...
) {
  .assertNamespace("ggplot2", "construct SSE plots")
  .assertSSEObject(x)
  method <- match.arg(method)
  nonpositivePolicy <- match.arg(nonpositivePolicy)

  if (identical(method, "distribution_mle")) {
    # Resolve comparisons EXACTLY ONCE here, unsuppressed, so the
    # inferred-df warning fires exactly once when it applies. It is tempting
    # to resolve again "for real" inside .ppePowerPlotData() below and
    # suppress this pass instead -- that was tried and is wrong: a second
    # resolve handed a non-NULL (already-resolved) comparisons list never
    # re-enters .resolveComparisons()'s comparisons-is-NULL branch, which is
    # the ONLY place .legacyComparisons() (and its warning) ever runs, so
    # the "real" pass would silently never warn at all rather than warning
    # twice. Splitting by mode up front also means a Type-I comparison never
    # reaches the sample-size curve and a power comparison never reaches the
    # point-range.
    mode_cmps <- .resolveComparisons(
      x, comparisons %||% x$runInfo$comparisons, models = models, ppe = TRUE
    )
    is_type1 <- vapply(mode_cmps, function(cmp) identical(cmp$mode, "type1"), logical(1))
    type1_cmps <- mode_cmps[is_type1]
    power_cmps <- mode_cmps[!is_type1]

    if (length(type1_cmps) > 0L && length(power_cmps) > 0L) {
      bad <- vapply(type1_cmps, `[[`, character(1), "label")
      cli::cli_warn(c(
        "!" = "{length(bad)} Type-I comparison{?s} omitted from the power curve: {.val {bad}}.",
        "i" = "A Type-I comparison estimates {.field df}, not a noncentrality that scales with study size, so it has no sample-size curve.",
        "i" = "See {.fn ppeSummary} for its estimated rate, or call {.fn plotSSEPpePower} with only Type-I comparisons for its point-range."
      ))
    }

    if (length(power_cmps) == 0L) {
      return(.plotSSEPpeType1PointRange(
        x, comparisons = type1_cmps, thresholds = thresholds,
        nonpositivePolicy = nonpositivePolicy, conf.level = conf.level,
        bootstrapSamples = bootstrapSamples, bootSeed = bootSeed
      ))
    }

    # power_cmps is threaded straight into .ppePowerPlotData() via
    # resolvedComparisons, NOT back through `comparisons` -- passing it
    # through `comparisons` would have .ppePowerPlotDataMle() resolve it a
    # second time, the exact mechanism that dropped the warning to zero
    # fires above. resolvedComparisons bypasses resolution there entirely.
    data <- .ppePowerPlotData(
      x,
      thresholds = thresholds,
      models = models,
      studySizes = studySizes,
      targetPower = targetPower,
      conf.level = conf.level,
      method = method,
      comparisons = NULL,
      nonpositivePolicy = nonpositivePolicy,
      bootstrapSamples = bootstrapSamples,
      bootSeed = bootSeed,
      resolvedComparisons = power_cmps
    )
  } else {
    data <- .ppePowerPlotData(
      x,
      thresholds = thresholds,
      models = models,
      studySizes = studySizes,
      targetPower = targetPower,
      conf.level = conf.level,
      method = method,
      comparisons = comparisons,
      nonpositivePolicy = nonpositivePolicy,
      bootstrapSamples = bootstrapSamples,
      bootSeed = bootSeed
    )
  }

  if (nrow(data) == 0L) {
    return(
      ggplot2::ggplot() +
        ggplot2::labs(title = "Parametric power estimation", subtitle = "No data to display") +
        ggplot2::theme_void()
    )
  }

  sample_info <- .studySampleSize(x)

  p <- ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = study_size,
      y = power,
      colour = model_label,
      fill = model_label,
      group = model_label
    )
  )

  # power_lower/power_upper are NA when bootstrapSamples = 0 (the caller
  # opted out) or when every bootstrap refit failed for a comparison; a
  # ribbon built from all-NA bounds is not a meaningful signal, so skip the
  # layer entirely rather than let ggplot2 draw nothing but warn about it.
  has_interval <- any(is.finite(data$power_lower) & is.finite(data$power_upper))
  if (has_interval) {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = power_lower, ymax = power_upper),
      alpha = 0.2,
      colour = NA
    )
  }

  p +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_vline(
      xintercept = sample_info$size,
      linetype = 2,
      linewidth = 0.5,
      colour = "grey40"
    ) +
    ggplot2::facet_wrap(
      ~ threshold,
      labeller = ggplot2::label_bquote(delta[OFV] == .(threshold))
    ) +
    ggplot2::labs(
      x = paste("Study size (", sample_info$unit, ")", sep = ""),
      y = "Estimated power (%)",
      colour = "Model",
      fill = "Model",
      title = "Parametric power estimation"
    ) +
    ggplot2::coord_cartesian(ylim = c(0, 100)) +
    ggplot2::theme_minimal()
}

#' Plot a multi-panel SSE diagnostics view
#'
#' @param x An `nlmixr2SSE` object.
#' @param biasStatistic Parameter-summary statistic to use in the bias panel.
#' @param ... Reserved for future extensions.
#'
#' @return A patchwork object.
#' @export
plotSSEDiagnostics <- function(
  x,
  biasStatistic = "relative_bias",
  ...
) {
  .assertNamespace("patchwork", "construct the diagnostics panel")
  .assertSSEObject(x)

  p1 <- plotSSEParameterBias(x, statistic = biasStatistic, ...)
  p2 <- plotSSEParameterEstimates(x, ...)
  p3 <- plotSSEOfvDistribution(x, ...)
  p4 <- plotSSEPower(x, ...)

  (p1 | p2) / (p3 | p4)
}
