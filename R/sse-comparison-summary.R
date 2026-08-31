# Paired empirical operating characteristics for explicit sseComparison()s.
#
# Two ideas drive this file:
#   1. Paired denominators. A replicate only counts if BOTH models produced a
#      finite, accepted OFV for that sample. Failed fits are excluded, never
#      imputed, and every count is reported so a reader can audit exclusions.
#   2. The test statistic comes from the comparison, not from a sign
#      convention: T = OFV_reduced - OFV_full, read off the comparison's own
#      members. That is what lets a Type-I comparison (simulation model is
#      the reduced one) come out positive.

# Per-model OFV vectors indexed by sample, filtered the same way the
# existing summary tables are (.computeOfvSummary(), .ofvDeltaPlotData()):
# non-accepted rows (error_message set, or excluded by outFilter) drop out,
# so a filtered-out replicate is treated as non-evaluable rather than
# missing-but-real.
.ofvByLabel <- function(x) {
  row_mask <- .summaryMask(x$rawResults, outFilter = x$runInfo$control$outFilter %||% NULL)
  filtered <- x$rawResults[row_mask, , drop = FALSE]

  # Anchor on the DECLARED replicate set, not on whatever reached rawResults.
  # A replicate where both models failed before writing any row would otherwise
  # vanish from n_attempted, n_paired_evaluable and n_excluded alike, so no
  # count would record that it ever existed -- the exact blind spot this
  # module's paired-denominator design exists to eliminate.
  observed <- sort(unique(x$rawResults$sample[x$rawResults$sample > 0L]))
  declared <- x$runInfo$samples
  samples <- if (!is.null(declared) && length(declared) == 1L && is.finite(declared)) {
    sort(union(seq_len(as.integer(declared)), observed))
  } else {
    observed
  }
  labels <- .knownModelLabels(x)
  ofv <- lapply(labels, function(label) {
    rows <- filtered$sample > 0L & filtered$model_label == label
    vec <- rep(NA_real_, length(samples))
    if (any(rows)) vec[match(filtered$sample[rows], samples)] <- filtered$objf[rows]
    vec
  })
  names(ofv) <- labels
  list(samples = samples, ofv = ofv)
}

.comparisonTestStatistic <- function(x, comparison) {
  by_label <- .ofvByLabel(x)
  full <- by_label$ofv[[comparison$full]]
  reduced <- by_label$ofv[[comparison$reduced]]
  # Computed from the comparison, never from a sign convention: this is what
  # makes a Type-I comparison (simulation model reduced) come out positive.
  stats_vec <- reduced - full
  stats_vec[is.finite(stats_vec)]
}

#' Paired-evaluable conditional rejection rate for explicit comparisons
#'
#' Computes a paired-evaluable rejection rate for one or more
#' [sseComparison()]s: for each comparison, the test statistic
#' `T = OFV_reduced - OFV_full` is computed per replicate, and a replicate
#' only counts toward the denominator if both models produced a finite,
#' accepted OFV for it. Uncertainty is an exact Clopper-Pearson binomial
#' interval on that empirical proportion -- purely a function of the observed
#' counts, with no distributional assumption about `T`. This is distinct from
#' the model-based (parametric bootstrap) uncertainty computed elsewhere,
#' hence `interval_type` and the `mcse_probability` naming.
#'
#' @details
#' The reported `probability` is the **paired-evaluable conditional
#' rejection probability**,
#' `P(T > criticalValue | both models produced a finite, accepted OFV)` --
#' not the unconditional rejection probability over every attempted
#' replicate. `n_paired_evaluable` (not `n_attempted`) is the denominator;
#' `n_excluded` counts replicates dropped because either fit failed or was
#' filtered out. If convergence or filtering depends on the simulated data,
#' the model, or the test statistic itself, that exclusion is informative
#' and the conditional rate can diverge from the unconditional one even when
#' the paired-evaluable fraction is high -- the `minPairedFraction` warning
#' only bounds how small that fraction is allowed to get, not how biased the
#' conditional rate can be at a given fraction. The unconditional design
#' operating characteristic is not identified without an explicit policy for
#' what an unpaired replicate should count as; when failures are
#' non-negligible, treat this rate as one input to a sensitivity analysis
#' under multiple failure policies, not as a final answer on its own.
#' Likewise `ci_lower`/`ci_upper` is an interval for the conditional
#' Bernoulli rate among retained pairs only -- it does not incorporate
#' uncertainty from the excluded replicates. Reserve "empirical power" and
#' "empirical Type I error" for contexts where this complete-case
#' conditioning has been checked and is acceptable.
#'
#' `mcse_probability` collapses to exactly 0 when `probability` is 0 or 1 --
#' precisely where the Clopper-Pearson interval is widest. At those
#' boundaries the interval (`ci_lower`/`ci_upper`), not `mcse_probability`,
#' is the honest summary of uncertainty; a reader who looks only at the
#' Monte Carlo standard error there would wrongly conclude the estimate is
#' exact.
#'
#' @param x An `nlmixr2SSE` object.
#' @param comparisons One or more [sseComparison()] objects, or `NULL` to use
#'   `x$runInfo$comparisons`.
#' @param models Optional subset of alternative-model labels, used only when
#'   `comparisons` is `NULL` to build the legacy simulation-vs-alternative
#'   comparisons. Mutually exclusive with `comparisons`.
#' @param conf.level Confidence level for the binomial interval.
#' @param minPairedFraction Minimum fraction of attempted replicates that must
#'   be paired evaluable before a low-yield warning is issued.
#' @return A `data.frame` with one row per comparison.
#' @export
comparisonSummary <- function(x, comparisons = NULL, models = NULL,
                              conf.level = 0.95, minPairedFraction = 0.5) {
  .assertSSEObject(x)
  cmps <- .resolveComparisons(x, comparisons %||% x$runInfo$comparisons, models)
  by_label <- .ofvByLabel(x)

  rows <- lapply(cmps, function(cmp) {
    full <- by_label$ofv[[cmp$full]]
    reduced <- by_label$ofv[[cmp$reduced]]
    paired <- is.finite(full) & is.finite(reduced)
    stat <- (reduced - full)[paired]
    n_paired <- sum(paired)
    if (n_paired == 0L) {
      .abortSSE("Comparison {.val {cmp$label}} has no paired evaluable replicates.")
    }
    n_exceeding <- sum(stat > cmp$criticalValue)
    interval <- .binomialInterval(n_exceeding, n_paired, conf.level = conf.level)
    prob <- n_exceeding / n_paired
    fraction <- n_paired / length(full)
    if (fraction < minPairedFraction) {
      # Raw counts, not just a rounded percentage: round(100 * 1/200) is 0, so a
      # percentage alone can report "0%" when replicates do remain, and "100%"
      # when some were excluded -- the opposite of what this warning is for.
      cli::cli_warn(c(
        "!" = "Only {n_paired} of {length(full)} replicates are paired evaluable for {.val {cmp$label}}.",
        "i" = "Failed fits are excluded, not imputed; operating characteristics may be biased."
      ))
    }
    data.frame(
      comparison = cmp$label, full = cmp$full, reduced = cmp$reduced,
      mode = cmp$mode, df = cmp$df %||% NA_real_, df_source = cmp$dfSource,
      alpha = cmp$alpha, critical_value = cmp$criticalValue,
      n_attempted = length(full),
      n_full_evaluable = sum(is.finite(full)),
      n_reduced_evaluable = sum(is.finite(reduced)),
      n_paired_evaluable = n_paired,
      n_excluded = length(full) - n_paired,
      n_exceeding = n_exceeding,
      probability = prob,
      mcse_probability = sqrt(prob * (1 - prob) / n_paired),
      ci_lower = interval[["lower"]], ci_upper = interval[["upper"]],
      interval_type = "empirical_binomial",
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}
