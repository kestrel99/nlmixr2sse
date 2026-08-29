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

.assertComparisonLabelsExist <- function(comparisons, labels, simLabel) {
  named <- unlist(lapply(comparisons, function(cmp) {
    resolve1 <- function(v) if (identical(v, "simulation")) simLabel else v
    c(resolve1(cmp$full), resolve1(cmp$reduced))
  }))
  unknown <- setdiff(unique(named), labels)
  if (length(unknown) > 0L) {
    .abortSSE(c(
      "{.arg comparisons} name{?s} model{?s} that this run will not fit: {.val {unknown}}.",
      "i" = "Models in this run: {.val {labels}}.",
      "i" = "Use {.val simulation} for the simulation model, or add the model to {.arg alternativeModels}."
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
  resolve1 <- function(value) if (identical(value, "simulation")) sim else value

  comparison$full <- resolve1(comparison$full)
  comparison$reduced <- resolve1(comparison$reduced)

  unknown <- setdiff(c(comparison$full, comparison$reduced), labels)
  if (length(unknown) > 0L) {
    .abortSSE(c(
      "Comparison {.val {comparison$label}} names unknown model{?s} {.val {unknown}}.",
      "i" = "Available labels: {.val {labels}}."
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
  labels <- vapply(out, `[[`, character(1), "label")
  if (anyDuplicated(labels) > 0L) {
    .abortSSE("Comparison labels must be unique; {.val {labels[duplicated(labels)]}} repeats.")
  }
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
