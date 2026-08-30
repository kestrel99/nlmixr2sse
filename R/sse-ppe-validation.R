#' Validate the proportional-noncentrality assumption across study sizes
#'
#' Parametric power estimation (PPE) extrapolates power to sample sizes that
#' were never simulated by assuming the noncentrality parameter scales
#' linearly with study size: `ncp = lambda_per_subject * subjects`. That
#' assumption holds for a simple, well-identified model but can fail when,
#' for example, the design or estimation method degrades disproportionately
#' at small sample sizes. This function checks the assumption empirically by
#' fitting the same comparison across multiple completed SSE runs that only
#' differ in study size, and reporting how consistent the fitted
#' noncentrality-per-subject ratio is across them.
#'
#' Only power comparisons (an estimated noncentrality `ncp` at a fixed `df`)
#' are supported: a Type-I comparison fixes `ncp = 0` and estimates `df`
#' instead, so there is no noncentrality for "proportional to study size" to
#' describe.
#'
#' @param runs A list of at least 2 completed `nlmixr2SSE` objects, each
#'   representing the same comparison run at a different study size.
#' @param comparisons Optional `sseComparison()` object (a single value,
#'   reused for every run) identifying which comparison to validate. When
#'   `NULL`, each run's own `runInfo$comparisons` is used (with a warning if
#'   `df` had to be inferred; see `.resolveComparisons()`).
#' @param conf.level Confidence level for each run's parametric-bootstrap
#'   interval.
#' @param bootstrapSamples Number of parametric-bootstrap replicates per run.
#' @param tolerance Fraction of the mean `lambda_per_subject` that the
#'   per-run values are allowed to spread over before the assumption is
#'   flagged as violated. Defaults to `0.25` (25%).
#' @return A list with:
#'   \describe{
#'     \item{table}{A `data.frame`, one row per run ordered by study size,
#'       with `subjects`, the fitted `lambda` (noncentrality estimate),
#'       its bootstrap `ci_lower`/`ci_upper`, and `lambda_per_subject`.}
#'     \item{lackOfFit}{The residual sum of squares from a through-origin fit
#'       of `lambda` on `subjects`, relative to the total sum of squares
#'       around the mean `lambda`. `NA` when fewer than 3 runs are supplied
#'       (a line always fits 2 points, so the statistic is uninformative).}
#'     \item{nonlinear}{`TRUE` when `lambda_per_subject` spreads over more
#'       than `tolerance` of its mean across runs.}
#'     \item{tolerance}{The tolerance used, echoed back for reference.}
#'   }
#' @export
validateSSEPpeScaling <- function(runs, comparisons = NULL, conf.level = 0.95,
                                  bootstrapSamples = 1000L, tolerance = 0.25) {
  if (!is.list(runs) || length(runs) < 2L) {
    .abortSSE("{.arg runs} must be a list of at least 2 completed SSE objects.")
  }
  lapply(runs, .assertSSEObject)

  fits <- lapply(runs, function(run) {
    resolved <- .resolveComparisons(run, comparisons %||% run$runInfo$comparisons,
                                    ppe = TRUE)
    if (length(resolved) != 1L) {
      .abortSSE(c(
        "Each run must resolve to exactly one comparison to validate.",
        "i" = "Got {length(resolved)} for one of the supplied runs.",
        "i" = "Pass a single {.fn sseComparison} via {.arg comparisons} to say which one to check."
      ))
    }
    cmp <- resolved[[1L]]
    if (!identical(cmp$mode, "power")) {
      .abortSSE(c(
        "{.arg comparisons} must be a power comparison (an estimated {.val ncp}).",
        "i" = "Comparison {.val {cmp$label}} is a Type-I comparison, which fixes {.val ncp = 0} and estimates {.val df} instead.",
        "i" = "Proportional noncentrality is a power-curve concept; there is no noncentrality to scale for a Type-I comparison."
      ))
    }
    subjects <- .studySampleSize(run)$size
    if (is.na(subjects) || subjects <= 0) {
      .abortSSE(c(
        "Every run needs a positive {.field runInfo$studySampleSize}.",
        "i" = "Found {.val {subjects}} for one of the supplied runs."
      ))
    }
    list(cmp = cmp,
         fit = .ppeFit(run, cmp, conf.level = conf.level,
                       bootstrapSamples = bootstrapSamples),
         subjects = subjects)
  })

  # Never pool replicates across study sizes, and never compare runs whose
  # definitions differ: a differing df or model label makes the noncentralities
  # incommensurable.
  keys <- vapply(fits, function(f) {
    paste(f$cmp$full, f$cmp$reduced, f$cmp$df, sep = "|")
  }, character(1))
  if (length(unique(keys)) > 1L) {
    .abortSSE(c(
      "The supplied runs have incompatible comparison definitions.",
      "i" = "Found: {.val {unique(keys)}}."
    ))
  }

  table <- data.frame(
    subjects = vapply(fits, `[[`, numeric(1), "subjects"),
    lambda = vapply(fits, function(f) f$fit$estimate, numeric(1)),
    ci_lower = vapply(fits, function(f) f$fit$ci_lower, numeric(1)),
    ci_upper = vapply(fits, function(f) f$fit$ci_upper, numeric(1)),
    stringsAsFactors = FALSE
  )
  table$lambda_per_subject <- table$lambda / table$subjects
  table <- table[order(table$subjects), , drop = FALSE]

  min_runs_for_lack_of_fit <- 3L
  lack_of_fit <- NA_real_
  if (nrow(table) >= min_runs_for_lack_of_fit) {
    # Through-origin line: lambda = beta * subjects.
    beta <- sum(table$lambda * table$subjects) / sum(table$subjects^2)
    resid <- table$lambda - beta * table$subjects
    lack_of_fit <- sum(resid^2) / sum((table$lambda - mean(table$lambda))^2)
  }
  spread <- diff(range(table$lambda_per_subject)) / mean(table$lambda_per_subject)
  nonlinear <- spread > tolerance
  if (nonlinear) {
    cli::cli_warn(c(
      "!" = "Estimated noncentrality is not proportional to study size across these runs.",
      "i" = "lambda per subject ranges over {round(100 * spread)}% of its mean.",
      "i" = "Extrapolating power from a single study size is unsupported over this range."
    ))
  }
  list(table = table, lackOfFit = lack_of_fit, nonlinear = nonlinear,
       tolerance = tolerance)
}
