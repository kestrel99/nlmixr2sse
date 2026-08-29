#' Define an explicit model comparison
#'
#' @param full,reduced Model labels, or the reserved token `"simulation"` which
#'   resolves to the label of the fitted simulation model.
#' @param df Degrees of freedom for the ordinary chi-square reference.
#' @param alpha Type-I error rate used to derive the critical value.
#' @param criticalValue An explicit critical value, for cases where the ordinary
#'   chi-square reference is inappropriate. Supplying it disables PPE.
#' @param label Optional unique comparison label.
#' @return An `sseComparison` object.
#' @export
sseComparison <- function(full, reduced, df = NULL, alpha = 0.05,
                          criticalValue = NULL, label = NULL) {
  checkmate::assertString(full, min.chars = 1L)
  checkmate::assertString(reduced, min.chars = 1L)
  if (identical(full, reduced)) {
    .abortSSE("{.arg full} and {.arg reduced} must be distinct model labels.")
  }
  if (is.null(df) == is.null(criticalValue)) {
    .abortSSE(c(
      "Supply exactly one of {.arg df} and {.arg criticalValue}.",
      "i" = "{.arg df} with {.arg alpha} gives the ordinary chi-square reference.",
      "i" = "{.arg criticalValue} is for references the chi-square does not cover."
    ))
  }
  checkmate::assertNumber(alpha, lower = 1e-10, upper = 0.5)
  if (!is.null(df)) {
    checkmate::assertNumber(df, lower = 1e-8, finite = TRUE)
    criticalValue <- stats::qchisq(1 - alpha, df = df)
  } else {
    checkmate::assertNumber(criticalValue, finite = TRUE)
  }
  structure(
    list(
      full = full,
      reduced = reduced,
      df = df,
      alpha = alpha,
      criticalValue = criticalValue,
      # PPE assumes a noncentral chi-square alternative, which a custom
      # critical value gives no basis for. Never let it be silently implied.
      ppeEligible = !is.null(df),
      dfSource = if (is.null(df)) NA_character_ else "explicit",
      label = label %||% sprintf("%s vs. %s", full, reduced)
    ),
    class = "sseComparison"
  )
}

.normalizeComparisonsArg <- function(comparisons) {
  if (is.null(comparisons)) {
    return(NULL)
  }
  if (inherits(comparisons, "sseComparison")) {
    comparisons <- list(comparisons)
  }
  lapply(comparisons, function(cmp) {
    if (!inherits(cmp, "sseComparison")) {
      .abortSSE("Each element of {.arg comparisons} must come from {.fn sseComparison}.")
    }
    cmp
  })
}

# One shared resolver. Duplicating this inline in `.resolveComparison()` and
# `.assertComparisonLabelsExist()` previously left the post-resolution
# distinctness check missing from one of the two call sites.
.resolveModelToken <- function(value, simLabel) {
  if (identical(value, "simulation")) simLabel else value
}

.assertComparisonLabelsUnique <- function(comparisons) {
  labels <- vapply(comparisons, `[[`, character(1), "label")
  if (anyDuplicated(labels) > 0L) {
    repeated <- unique(labels[duplicated(labels)])
    .abortSSE("Comparison labels must be unique; {.val {repeated}} {?repeats/repeat}.")
  }
}

.assertComparisonLabelsExist <- function(comparisons, labels, simLabel) {
  resolved <- lapply(comparisons, function(cmp) {
    list(
      label = cmp$label,
      full = .resolveModelToken(cmp$full, simLabel),
      reduced = .resolveModelToken(cmp$reduced, simLabel)
    )
  })

  named <- unlist(lapply(resolved, function(r) c(r$full, r$reduced)))
  unknown <- setdiff(unique(named), labels)
  if (length(unknown) > 0L) {
    .abortSSE(c(
      "Unknown model{?s} in {.arg comparisons}: {.val {unknown}}.",
      "i" = "Models in this run: {.val {labels}}.",
      "i" = "Use {.val simulation} for the simulation model, or add the model to {.arg alternativeModels}."
    ))
  }

  # The constructor compares raw strings, so it cannot catch
  # sseComparison(full = "<sim label>", reduced = "simulation"): both sides
  # resolve to the same model. Left unchecked that yields a model compared
  # against itself, every test statistic identically zero, and a meaningless
  # PPE fit reported as if it were real.
  self_compared <- Filter(function(r) identical(r$full, r$reduced), resolved)
  if (length(self_compared) > 0L) {
    bad <- self_compared[[1L]]
    .abortSSE(c(
      "Comparison {.val {bad$label}} resolves both members to {.val {bad$full}}.",
      "i" = "A model compared against itself has no degrees of freedom and no test statistic.",
      "i" = "The {.val simulation} token resolves to {.val {simLabel}}."
    ))
  }
}

.simulationLabel <- function(x) .simulationSpec(x)$label

.knownModelLabels <- function(x) {
  vapply(x$runInfo$fitSpecs %||% list(), `[[`, character(1), "label")
}

.resolveComparison <- function(x, comparison) {
  labels <- .knownModelLabels(x)
  sim <- .simulationLabel(x)

  comparison$full <- .resolveModelToken(comparison$full, sim)
  comparison$reduced <- .resolveModelToken(comparison$reduced, sim)

  unknown <- setdiff(c(comparison$full, comparison$reduced), labels)
  if (length(unknown) > 0L) {
    .abortSSE(c(
      "Comparison {.val {comparison$label}} names unknown model{?s} {.val {unknown}}.",
      "i" = "Available labels: {.val {labels}}."
    ))
  }

  # The constructor compares raw strings, so it cannot catch
  # sseComparison(full = "<sim label>", reduced = "simulation"): both sides
  # resolve to the same model. Left unchecked that yields a model compared
  # against itself, every test statistic identically zero, and a meaningless
  # PPE fit reported as if it were real.
  if (identical(comparison$full, comparison$reduced)) {
    .abortSSE(c(
      "Comparison {.val {comparison$label}} resolves both members to {.val {comparison$full}}.",
      "i" = "A model compared against itself has no degrees of freedom and no test statistic.",
      "i" = "The {.val simulation} token resolves to {.val {sim}}."
    ))
  }

  # The mode follows from which member was simulated: that is the model whose
  # hypothesis is true, so it decides whether the run measures power or Type-I.
  comparison$mode <- if (identical(comparison$full, sim)) {
    "power"
  } else if (identical(comparison$reduced, sim)) {
    "type1"
  } else {
    .abortSSE(c(
      "Comparison {.val {comparison$label}} names neither member as the simulation model.",
      "i" = "The simulation model is {.val {sim}}.",
      "i" = "Without it, no hypothesis is known true and neither power nor Type-I is defined."
    ))
  }
  comparison
}

.resolveComparisons <- function(x, comparisons = NULL, models = NULL, ppe = FALSE) {
  if (!is.null(comparisons) && !is.null(models)) {
    .abortSSE(c(
      "Supply either {.arg comparisons} or {.arg models}, not both.",
      "i" = "{.arg comparisons} already names every pair explicitly."
    ))
  }
  if (is.null(comparisons)) {
    comparisons <- .legacyComparisons(x, models = models, ppe = ppe)
  }
  if (inherits(comparisons, "sseComparison")) comparisons <- list(comparisons)
  out <- lapply(comparisons, function(cmp) .resolveComparison(x, cmp))
  .assertComparisonLabelsUnique(out)
  out
}

.legacyComparisons <- function(x, models = NULL, ppe = FALSE) {
  sim <- .simulationLabel(x)
  specs <- Filter(function(s) identical(s$role, "alternative"), x$runInfo$fitSpecs %||% list())
  alt_labels <- vapply(specs, `[[`, character(1), "label")
  if (!is.null(models)) alt_labels <- intersect(alt_labels, models)
  if (ppe) {
    cli::cli_warn(c(
      "!" = "Degrees of freedom inferred from parameter counts for {length(alt_labels)} comparison{?s}.",
      "i" = "Parametric power estimation treats df as known; an inferred value is a convenience, not an assertion.",
      "i" = "Define comparisons explicitly with {.fn sseComparison} to remove this warning."
    ))
  }
  lapply(alt_labels, function(label) {
    df <- .modelDegreesFreedom(x, label)
    cmp <- sseComparison(sim, label, df = df)
    cmp$dfSource <- "parameter_count"
    cmp
  })
}
