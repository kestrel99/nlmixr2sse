# SSE Statistical Rigour Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Strengthen the statistical definition, diagnostics, and uncertainty
quantification of `nlmixr2sse` — explicit comparisons, paired denominators,
role-specific estimation starts, distribution-based PPE with bootstrap
uncertainty and adequacy diagnostics, and parameter-draw adequacy summaries.

**Architecture:** Preserve the existing simulation and estimation engine and the
`independent_iw` and `joint` parameter-draw modes. Add an explicit comparison
object that is defined at run time and persisted; compute every model-comparison
result on paired evaluable replicates; replace threshold-by-threshold PPE as the
default with a full-distribution noncentral chi-square likelihood fit plus
parametric-bootstrap uncertainty, distributional diagnostics, and sample-size
scaling validation. Retain current behaviour behind explicit legacy options.

**Tech Stack:** R package; `nlmixr2est`, `rxode2`, `nlmixr2utils` (>= 0.3),
`ggplot2`, `checkmate`, `cli`, `lifecycle`, testthat 3e, base R `stats`.

---

## Source documents and reconciliation

This plan merges two documents:

- `docs/superpowers/plans/2026-08-29-sse-statistical-rigour.md` — the parent
  roadmap, 11 tasks. Its task order and numbering are preserved here.
- `docs/superpowers/specs/2026-08-30-sse-psn-ppe-design.md` — the detailed PPE
  design, folded into Tasks 2 and 5–8.

Where they disagreed, the resolution is:

| Question | Resolution |
| --- | --- |
| Default PPE method | `distribution_mle`. **Breaking change to plot output**; needs a NEWS entry. |
| Method names | `distribution_mle` / `exceedance`. Not `psn_mle`: this package does not seek PsN execution, seed, or file-format parity, only the same estimator. |
| `sseComparison()` signature | The rigour superset: `(full, reduced, df, alpha, criticalValue, label)`, definable at run time and persisted in `runInfo`. |
| Diagnostics function | `plotSSEPpeDiagnostics()`, plural, matching the existing `plotSSEDiagnostics()`. |
| Bootstrap argument | `bootstrapSamples`, not `nBoot`. |
| Non-positive ΔOFV policy | The rigour `nonpositive = c("warn", "error", "drop")`; counts are never suppressed. |
| Summary functions | Both. `comparisonSummary()` is empirical (counts, binomial CI); `ppeSummary()` is model-based (MLE, bootstrap CI). Rigour decision 9 requires distinct names. |
| Threshold handling | The PPE spec kept per-threshold faceting under the MLE. Superseded: under `distribution_mle` the critical value comes from the comparison (`qchisq(1 - alpha, df)`, or an explicit `criticalValue`), so one comparison yields one curve. Faceting over a `thresholds` vector was an artefact of the exceedance estimator needing a threshold, and is retained only for `method = "exceedance"`. |

### Findings verified before planning

These were reproduced directly and drive several design choices:

- **Truncation.** The estimator drops ΔOFV ≤ 0 and fits the *unconditional*
  noncentral chi-square density. Since `P(X > 0) = 1` for that distribution, the
  discarded values are ones the model says cannot occur, and dropping them
  biases ncp upward. Counts must always be reported (Task 5).
- **Boundary solutions are common.** When the retained mean falls below `df`,
  the constrained MLE sits at the lower bound and estimated power equals
  `alpha` exactly. Measured: 0/200 replicates at ncp=1/df=5/n=50, 1/200 at
  ncp=2/df=8/n=40, and **31/200** at ncp=0.5/df=4/n=60. Must be flagged, not
  silently plotted (Task 5).
- **An out-of-bounds `optim()` start is harmless.** `optim()` projects it into
  the feasible region and converges normally (`init = -3` on `(p-5)^2` returns
  `par = 5`, `convergence = 0`). No guard is warranted.
- **The parametric bootstrap does not reproduce the truncation.** `rchisq()`
  draws are strictly positive, so no replicate ever discards a value while the
  real data did. The interval is model-based only, and must be named as such
  (Task 6).
- **Sign-convention trap.** `.ofvDeltaPlotData()` computes
  `OFV(simulation) - OFV(alternative)` (`R/plot-sse.R:203`) and
  `.ppePowerPlotData()` negates it, so the statistic is
  `OFV(alternative) - OFV(simulation)`. That is non-negative only when the
  simulation model is the full one; under a Type-I run every value goes
  negative and truncation discards all of them. Tasks 2 and 3 compute
  `OFV_reduced - OFV_full` from the comparison instead.
- **df inference is silently fragile.** `.modelDegreesFreedom()`
  (`R/plot-sse.R:265`) returns `1L` when a label is missing from `fitSpecs` and
  clamps with `max(df, 1L)`. Neither path emits a signal (Task 2).

## File Structure

| File | Responsibility |
| --- | --- |
| `R/sse-comparison.R` *(new)* | `sseComparison()`, validation, role-token resolution, mode derivation, df precedence. |
| `R/sse-comparison-summary.R` *(new)* | Paired evaluable joins, empirical operating characteristics, binomial intervals. |
| `R/sse-ppe.R` *(new)* | Noncentral chi-square likelihood, MLE, parametric bootstrap, `ppeSummary()`. No `ggplot2`. |
| `R/plot-sse-ppe-diagnostics.R` *(new)* | `plotSSEPpeDiagnostics()`: ECDF, QQ, envelopes, discrepancy statistic. |
| `R/sse-ppe-validation.R` *(new)* | `validateSSEPpeScaling()` across study sizes. |
| `R/sse-parameter-diagnostics.R` *(new)* | `parameterDrawSummary()`, `plotSSEParameterDraws()`. |
| `R/plot-sse.R` | PPE curve rendering only; estimation moves out. |
| `R/sse-control.R` | `referenceInitials`, `alternativeInitials`, `randomEstimationInits` deprecation. |
| `R/run-sse.R`, `R/recompute-sse.R` | Accept, persist, and validate comparisons and start policy. |
| `R/sse-helpers.R`, `R/sse-methods.R` | MCSE fields, summary plumbing. |
| `tests/testthat/helper-ppe-fixtures.R` *(new)* | Deterministic ΔOFV and SSE fixtures for PPE work. |

Estimation is deliberately separated from rendering: the statistics are the
risky part and must be testable without constructing a plot.

---

## Task 1: Freeze the current statistical contract

**Files:**
- Create: `tests/testthat/test-sse-statistical-contract.R`
- Create: `tests/testthat/helper-ppe-fixtures.R`
- Read: `R/plot-sse.R:151-350`, `R/sse-helpers.R`

- [ ] **Step 1: Record the baseline test result**

Run: `NOT_CRAN=true Rscript -e 'pkgload::load_all("."); testthat::test_dir("tests/testthat", reporter="progress", stop_on_failure=FALSE)'`

Expected: `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 537 ]`. Record this number in the
commit message. Any later task that changes it must say why.

- [ ] **Step 2: Write the shared fixtures**

Every later task builds on these. The existing `fake_sse_object()` has only two
samples and fixed OFVs, which is far too little to fit a distribution to, so
these take the OFVs as arguments.

```r
# tests/testthat/helper-ppe-fixtures.R

# Deterministic noncentral chi-square draws that leave the caller's RNG alone.
ppe_dofv <- function(n = 200L, df = 1, ncp = 8, seed = 101L) {
  has <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old <- if (has) get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit({
    if (has) {
      assign(".Random.seed", old, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      # set.seed() creates .Random.seed where there was none; restoring the
      # caller's state means restoring its absence too.
      rm(".Random.seed", envir = globalenv())
    }
  }, add = TRUE)
  set.seed(seed)
  stats::rchisq(n, df = df, ncp = ncp)
}

ppe_initial_estimates <- function(n) {
  data.frame(
    sample = seq_len(n),
    tka = seq(0.50, 0.70, length.out = n),
    tcl = seq(1.00, 1.20, length.out = n),
    "omega(eta.ka,eta.ka)" = seq(0.30, 0.32, length.out = n),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

# An nlmixr2SSE whose per-sample OFVs are exactly the supplied vectors. NA
# entries become missing fits, so asymmetric evaluability is exercisable.
fake_paired_sse_object <- function(full_ofv, reduced_ofv,
                                   fit_label = "fake_sse_fit",
                                   alt_label = "alt1",
                                   subjects = 12L, seed = 42L) {
  stopifnot(length(full_ofv) == length(reduced_ofv))
  # Fail here rather than several frames down in .computeSSEOutputs(), which
  # reports only "argument is of length zero" when handed no rows at all.
  stopifnot(
    "need at least one finite OFV" =
      any(is.finite(full_ofv) | is.finite(reduced_ofv))
  )
  fit <- fake_sse_fit()
  n <- length(full_ofv)

  rows <- list()
  for (i in seq_len(n)) {
    if (is.finite(full_ofv[[i]])) {
      rows[[length(rows) + 1L]] <- nlmixr2utils::rawResultsRow(
        fit, source = "sse", hypothesis = "simulation", sample = i,
        modelLabel = fit_label, role = "simulation",
        theta = c(tka = 0.55, tcl = 1.05),
        omega = c("omega(eta.ka,eta.ka)" = 0.31),
        objf = full_ofv[[i]]
      )
    }
    if (is.finite(reduced_ofv[[i]])) {
      rows[[length(rows) + 1L]] <- nlmixr2utils::rawResultsRow(
        fit, source = "sse", hypothesis = "alternative_1", sample = i,
        modelLabel = alt_label, role = "alternative",
        theta = c(tka = 0.52, tcl = 1.04),
        objf = reduced_ofv[[i]]
      )
    }
  }
  raw_results <- do.call(rbind, rows)

  run_info <- list(
    # The real seed must be carried: Task 6's .ppeDefaultSeed() derives the
    # bootstrap seed from runInfo$seed, so a hardcoded value would make every
    # fixture bootstrap identically regardless of the data it was built from.
    fitName = fit_label, samples = n, seed = as.integer(seed),
    parameterSource = "fixed", estimateSimulation = TRUE,
    studySampleSize = subjects, studySampleUnit = "subjects",
    studyIdColumn = "ID", studyObservationCount = 2L * subjects,
    fitSpecs = fake_sse_fit_specs(fit_label, alt_label),
    control = runSSEControl(workers = 1L)
  )
  outputs <- .computeSSEOutputs(
    rawResults = raw_results,
    initialWide = ppe_initial_estimates(n),
    fitSpecsSnapshot = run_info$fitSpecs,
    runInfo = run_info
  )
  .newNlmixr2SSE(
    runInfo = run_info, rawResults = raw_results,
    alternativeSpecs = list(list(
      label = alt_label, est = "focei", control = list(print = 0L),
      isFit = TRUE, hasDataOverride = FALSE
    )),
    outputDir = tempdir(), timestamp = Sys.time(),
    referenceValues = outputs$referenceValues,
    initialValues = outputs$initialValues,
    parameterSummary = outputs$parameterSummary,
    ofvSummary = outputs$ofvSummary,
    powerSummary = outputs$powerSummary
  )
}

# An SSE object whose test statistics are exactly ppe_dofv() draws. Holding the
# full-model OFV constant makes OFV_reduced - OFV_full equal the draw exactly,
# so the estimator can be checked against a known noncentrality parameter.
fake_ppe_sse_object <- function(df = 1, ncp = 8, n = 200L, seed = 101L,
                                subjects = 12L, nNonPositive = 0L) {
  stopifnot(nNonPositive <= n)
  d <- ppe_dofv(n = n, df = df, ncp = ncp, seed = seed)
  if (nNonPositive > 0L) {
    idx <- seq_len(nNonPositive)
    d[idx] <- -abs(d[idx])
  }
  full <- rep(1000, n)
  fake_paired_sse_object(full_ofv = full, reduced_ofv = full + d,
                         subjects = subjects, seed = seed)
}
```

- [ ] **Step 3: Verify the fixtures build and behave**

```r
test_that("the paired fixture reproduces the supplied OFVs exactly", {
  sse <- fake_paired_sse_object(full_ofv = c(100, 101), reduced_ofv = c(105, 108))
  expect_equal(-.ofvDeltaPlotData(sse)$delta_ofv, c(5, 7))
})

test_that("the PPE fixture reproduces its draws as test statistics", {
  d <- ppe_dofv(n = 20L, df = 1, ncp = 8, seed = 5L)
  sse <- fake_ppe_sse_object(df = 1, ncp = 8, n = 20L, seed = 5L)
  expect_equal(-.ofvDeltaPlotData(sse)$delta_ofv, d)
})

test_that("ppe_dofv leaves the caller's RNG untouched", {
  set.seed(1L)
  before <- .Random.seed
  ppe_dofv(n = 10L, seed = 2L)
  expect_equal(.Random.seed, before)
})
```

Run: `NOT_CRAN=true Rscript -e 'pkgload::load_all("."); testthat::test_file("tests/testthat/test-sse-statistical-contract.R")'`
Expected: PASS.

- [ ] **Step 4: Write characterisation tests for the current behaviour**

```r
test_that("current delta sign convention is reference minus alternative", {
  sse <- fake_sse_object()
  delta <- .ofvDeltaPlotData(sse)

  # simulation objf 100/101, alternative objf 103/105
  expect_equal(delta$delta_ofv, c(-3, -4))
  # and the PPE test statistic is the negation
  expect_equal(-delta$delta_ofv, c(3, 4))
})

test_that("current PPE inverts the exceedance probability per threshold", {
  sse <- fake_sse_object()
  test_stat <- c(3, 4)  # -delta_ofv for the two-sample fixture

  low <- .ppePowerPlotData(sse, thresholds = 1, studySizes = 12L)
  high <- .ppePowerPlotData(sse, thresholds = 3.5, studySizes = 12L)

  # At the base study size the scaling factor is 1, so each threshold's ncp is
  # solved to reproduce that threshold's own clipped exceedance rate exactly.
  # A single ncp fitted to the whole distribution cannot reproduce two
  # different rates at once, so these equalities are precisely what breaks when
  # the estimator changes in Task 5. Asserting only that the two powers differ
  # would be useless: power depends on the threshold under any implementation.
  expect_equal(low$power, 100 * .clipProbability(mean(test_stat > 1), 2))
  expect_equal(high$power, 100 * .clipProbability(mean(test_stat > 3.5), 2))
  expect_equal(unique(low$degrees_freedom), 1L)
})

test_that("modelDegreesFreedom falls back to 1 for an unknown label", {
  sse <- fake_sse_object()
  expect_equal(.modelDegreesFreedom(sse, "no_such_model"), 1L)
})
```

- [ ] **Step 5: Run the characterisation tests**

Run: `NOT_CRAN=true Rscript -e 'pkgload::load_all("."); testthat::test_file("tests/testthat/test-sse-statistical-contract.R")'`
Expected: PASS. These describe today's behaviour, so they pass now and fail only
when a later task deliberately changes it.

- [ ] **Step 6: Commit**

```bash
git add tests/testthat/test-sse-statistical-contract.R tests/testthat/helper-ppe-fixtures.R
git commit -m "test: freeze the current SSE statistical contract and add PPE fixtures"
```

---

## Task 2: Explicit comparison objects

**Files:**
- Create: `R/sse-comparison.R`
- Create: `tests/testthat/test-sse-comparison.R`
- Modify: `R/run-sse.R`, `R/recompute-sse.R`, `NAMESPACE`

- [ ] **Step 1: Write the failing constructor tests**

```r
test_that("sseComparison requires distinct labels and one reference definition", {
  expect_error(sseComparison("a", "a", df = 1), "distinct")
  expect_error(sseComparison("a", "b"), "exactly one")
  expect_error(sseComparison("a", "b", df = 1, criticalValue = 3.84), "exactly one")
})

test_that("sseComparison derives the chi-square critical value", {
  cmp <- sseComparison("simulation", "no_cov", df = 1, alpha = 0.05)

  expect_s3_class(cmp, "sseComparison")
  expect_equal(cmp$criticalValue, stats::qchisq(0.95, df = 1))
  expect_equal(cmp$dfSource, "explicit")
  expect_equal(cmp$label, "simulation vs. no_cov")
})

test_that("an explicit criticalValue does not assert PPE validity", {
  cmp <- sseComparison("a", "b", criticalValue = 5)

  expect_equal(cmp$criticalValue, 5)
  expect_null(cmp$df)
  expect_false(cmp$ppeEligible)
})
```

- [ ] **Step 2: Run to verify failure**

Run: `NOT_CRAN=true Rscript -e 'pkgload::load_all("."); testthat::test_file("tests/testthat/test-sse-comparison.R")'`
Expected: FAIL, `could not find function "sseComparison"`.

- [ ] **Step 3: Implement the constructor**

```r
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
```

- [ ] **Step 4: Run to verify the constructor tests pass**

Run: `NOT_CRAN=true Rscript -e 'pkgload::load_all("."); testthat::test_file("tests/testthat/test-sse-comparison.R")'`
Expected: PASS.

- [ ] **Step 5: Write failing tests for resolution and mode derivation**

```r
test_that("the simulation role token resolves to the fitted label", {
  sse <- fake_sse_object()
  cmp <- .resolveComparison(sse, sseComparison("simulation", "alt1", df = 1))

  expect_equal(cmp$full, "fake_sse_fit")
  expect_equal(cmp$mode, "power")
})

test_that("a reduced simulation model gives type1 mode", {
  sse <- fake_sse_object()
  cmp <- .resolveComparison(sse, sseComparison("alt1", "simulation", df = 1))

  expect_equal(cmp$mode, "type1")
})

test_that("a comparison with neither member simulated is rejected", {
  sse <- fake_sse_object(alt_label = "alt1")
  expect_error(
    .resolveComparison(sse, sseComparison("alt1", "alt1x", df = 1)),
    "neither"
  )
})

test_that("an unknown label lists the available labels", {
  sse <- fake_sse_object()
  expect_error(
    .resolveComparison(sse, sseComparison("simulation", "nope", df = 1)),
    "fake_sse_fit"
  )
})

test_that("legacy comparisons mark df as inferred and warn once", {
  sse <- fake_sse_object()
  expect_warning(
    cmps <- .resolveComparisons(sse, comparisons = NULL, ppe = TRUE),
    "inferred from parameter counts"
  )
  expect_equal(cmps[[1L]]$dfSource, "parameter_count")
  expect_equal(cmps[[1L]]$full, "fake_sse_fit")
  expect_equal(cmps[[1L]]$reduced, "alt1")
})
```

- [ ] **Step 6: Run to verify failure**

Expected: FAIL, `could not find function ".resolveComparison"`.

- [ ] **Step 7: Implement resolution and mode derivation**

```r
.simulationLabel <- function(x) .simulationSpec(x)$label

.knownModelLabels <- function(x) {
  vapply(x$runInfo$fitSpecs %||% list(), `[[`, character(1), "label")
}

# One shared resolver. Duplicating this inline is what previously left the
# post-resolution distinctness check missing from one of two call sites.
.resolveModelToken <- function(value, simLabel) {
  if (identical(value, "simulation")) simLabel else value
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
  labels <- vapply(out, `[[`, character(1), "label")
  if (anyDuplicated(labels) > 0L) {
    repeated <- unique(labels[duplicated(labels)])
    .abortSSE("Comparison labels must be unique; {.val {repeated}} {?repeats/repeat}.")
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
```

- [ ] **Step 8: Run to verify the resolution tests pass**

Expected: PASS.

- [ ] **Step 9: Accept and persist comparisons in `runSSE()`**

In `R/run-sse.R`, add a `comparisons = NULL` argument, validate each element is
an `sseComparison`, and persist normalised copies:

```r
  run_info$comparisons <- if (is.null(comparisons)) {
    # Preserve, never wipe: a comparison is a reporting definition, so a resume
    # or addModels call that omits the argument must keep what was recorded,
    # exactly as parameterSourceInfo and the studySample* fields do.
    existing_run_info$comparisons
  } else {
    if (inherits(comparisons, "sseComparison")) comparisons <- list(comparisons)
    lapply(comparisons, function(cmp) {
      if (!inherits(cmp, "sseComparison")) {
        .abortSSE("Each element of {.arg comparisons} must come from {.fn sseComparison}.")
      }
      cmp
    })
  }
```

Validate the labels once the full fit-spec list is assembled, aborting on any
that will not exist. Every label is knowable at this point — the simulation
model plus every alternative, including those being added by `addModels` — so
deferring the check to reporting time would let a typo persist silently through
an entire run, which is the failure class this plan exists to remove. The
reserved `"simulation"` token is resolved before checking.

```r
  .assertComparisonLabelsExist <- function(comparisons, labels, simLabel) {
    named <- unlist(lapply(comparisons, function(cmp) {
      resolve1 <- function(v) if (identical(v, "simulation")) simLabel else v
      c(resolve1(cmp$full), resolve1(cmp$reduced))
    }))
    unknown <- setdiff(unique(named), labels)
    if (length(unknown) > 0L) {
      # Phrased to avoid subject-verb agreement entirely: cli's {?s} appends as
      # the quantity grows, which is right for a noun but inverted for a verb
      # ("name" for one, "names" for many is backwards).
      .abortSSE(c(
        "Unknown model{?s} in {.arg comparisons}: {.val {unknown}}.",
        "i" = "Models in this run: {.val {labels}}.",
        "i" = "Use {.val simulation} for the simulation model, or add the model to {.arg alternativeModels}."
      ))
    }
  }
```

This check must run on **every** path that accepts `comparisons`, including the
early return for an already-completed run directory. That path returns before
the main body executes, so validating only in the main body silently discards
the argument and lets a typo through with no error at all — the opposite of the
documented contract.

In `R/recompute-sse.R`, allow `comparisons` to replace `run_info$comparisons`
without refitting, since a comparison is a reporting definition rather than an
input to simulation or estimation.

- [ ] **Step 10: Run the full suite and commit**

Run: `NOT_CRAN=true Rscript -e 'pkgload::load_all("."); testthat::test_dir("tests/testthat", reporter="progress", stop_on_failure=FALSE)'`
Expected: `FAIL 0`, PASS count above 537.

```bash
git add R/sse-comparison.R tests/testthat/test-sse-comparison.R R/run-sse.R R/recompute-sse.R NAMESPACE
git commit -m "feat: explicit sseComparison objects with derived power/type1 mode"
```

---

## Task 3: Paired empirical operating characteristics

**Files:**
- Create: `R/sse-comparison-summary.R`
- Create: `tests/testthat/test-sse-comparison-summary.R`
- Modify: `R/sse-methods.R`

- [ ] **Step 1: Write the failing denominator test**

```r
test_that("the denominator is paired evaluable replicates, not any single model", {
  sse <- fake_paired_sse_object(
    full_ofv    = c(100, 101, 102, NA,  104),
    reduced_ofv = c(110, 111, NA,  113, 114)
  )
  cmp <- sseComparison("simulation", "alt1", df = 1)

  s <- comparisonSummary(sse, comparisons = cmp)

  expect_equal(s$n_attempted, 5L)
  expect_equal(s$n_full_evaluable, 4L)
  expect_equal(s$n_reduced_evaluable, 4L)
  expect_equal(s$n_paired_evaluable, 3L)   # samples 1, 2, 5
  expect_equal(s$n_excluded, 2L)
})

test_that("the empirical interval is reproducible from the reported counts", {
  sse <- fake_paired_sse_object(
    full_ofv    = rep(100, 20),
    reduced_ofv = 100 + c(rep(10, 12), rep(0.5, 8))
  )
  cmp <- sseComparison("simulation", "alt1", df = 1)

  s <- comparisonSummary(sse, comparisons = cmp)

  expect_equal(s$n_exceeding, 12L)
  expect_equal(s$probability, 12 / 20)
  expect_equal(
    c(s$ci_lower, s$ci_upper),
    c(stats::qbeta(0.025, 12, 9), stats::qbeta(0.975, 13, 8))
  )
})

test_that("the test statistic is reduced minus full regardless of which was simulated", {
  sse <- fake_paired_sse_object(full_ofv = c(100, 100), reduced_ofv = c(105, 108))

  power_cmp <- sseComparison("simulation", "alt1", df = 1)
  expect_equal(.comparisonTestStatistic(sse, .resolveComparison(sse, power_cmp)), c(5, 8))
})
```

- [ ] **Step 2: Run to verify failure**

Expected: FAIL, `could not find function "comparisonSummary"`.

- [ ] **Step 3: Implement the paired join and summary**

```r
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

#' Empirical operating characteristics for explicit comparisons
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
```

- [ ] **Step 4: Run to verify the tests pass**

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add R/sse-comparison-summary.R tests/testthat/test-sse-comparison-summary.R R/sse-methods.R
git commit -m "feat: paired-evaluable comparison summaries with binomial intervals"
```

---

## Task 4: Role-specific estimation starts

**Files:**
- Modify: `R/sse-control.R`, `R/sse-helpers.R`, `R/run-sse.R`, `R/recompute-sse.R`, `DESCRIPTION`
- Modify: `tests/testthat/test-control.R`, `tests/testthat/test-run-sse.R`

- [ ] **Step 1: Write the failing control tests**

```r
test_that("initials are settable per model role", {
  ctl <- runSSEControl(referenceInitials = "simulation", alternativeInitials = "model")

  expect_equal(ctl$referenceInitials, "simulation")
  expect_equal(ctl$alternativeInitials, "model")
})

test_that("initials reject unknown policies", {
  expect_error(runSSEControl(referenceInitials = "truth"), "must be one of")
})

test_that("randomEstimationInits maps to both roles with a deprecation warning", {
  expect_warning(
    ctl <- runSSEControl(randomEstimationInits = TRUE),
    "deprecated"
  )
  expect_equal(ctl$referenceInitials, "simulation")
  expect_equal(ctl$alternativeInitials, "simulation")
})

test_that("the old and new arguments cannot contradict each other", {
  expect_error(
    suppressWarnings(runSSEControl(randomEstimationInits = TRUE,
                                   referenceInitials = "model")),
    "contradict"
  )
})
```

- [ ] **Step 2: Run to verify failure**

Expected: FAIL, `unused argument (referenceInitials = "simulation")`.

- [ ] **Step 3: Add `lifecycle` and implement the arguments**

```bash
Rscript -e 'usethis::use_lifecycle()'
```

In `R/sse-control.R`:

```r
  referenceInitials = c("model", "simulation"),
  alternativeInitials = c("model", "simulation"),
  randomEstimationInits = lifecycle::deprecated(),
```

```r
  referenceInitials <- match.arg(referenceInitials)
  alternativeInitials <- match.arg(alternativeInitials)
  if (lifecycle::is_present(randomEstimationInits)) {
    lifecycle::deprecate_soft(
      "0.4.0", "runSSEControl(randomEstimationInits)",
      details = "Use `referenceInitials` and `alternativeInitials` instead."
    )
    mapped <- if (isTRUE(randomEstimationInits)) "simulation" else "model"
    if (!referenceInitialsMissing && !identical(referenceInitials, mapped)) {
      .abortSSE(
        "{.arg randomEstimationInits} and {.arg referenceInitials} contradict each other."
      )
    }
    if (!alternativeInitialsMissing && !identical(alternativeInitials, mapped)) {
      .abortSSE(
        "{.arg randomEstimationInits} and {.arg alternativeInitials} contradict each other."
      )
    }
    referenceInitials <- mapped
    alternativeInitials <- mapped
  }
```

`referenceInitialsMissing` and `alternativeInitialsMissing` are captured with
`missing()` at the top of the function body, before `match.arg()` runs.

- [ ] **Step 4: Apply the policy per role in `.fitTaskRecord()`**

```r
  initials_policy <- if (identical(spec$role, "simulation")) {
    control$referenceInitials
  } else {
    control$alternativeInitials
  }
  use_simulation_starts <- identical(initials_policy, "simulation")
```

Replace the existing `isTRUE(control$randomEstimationInits)` branch with
`use_simulation_starts`.

- [ ] **Step 5: Write the invariance test**

```r
test_that("starting-value policy changes no generating value or dataset", {
  base <- runSSE(fit, samples = 2L, seed = 7,
                 control = runSSEControl(workers = 1L, referenceInitials = "model"),
                 outputDir = dir_a, restart = TRUE)
  alt <- runSSE(fit, samples = 2L, seed = 7,
                control = runSSEControl(workers = 1L, referenceInitials = "simulation"),
                outputDir = dir_b, restart = TRUE)

  expect_equal(base$initialValues, alt$initialValues)
})
```

- [ ] **Step 6: Record and validate on resume**

Add `run_info$referenceInitials` and `run_info$alternativeInitials`, and extend
`.validateResumeRequest()` to abort on a mismatch, mirroring the existing
`rxThreads` check. Legacy directories recording only `randomEstimationInits`
resolve unambiguously to both roles, so unlike an unknown historical thread
count they are mapped rather than warned about.

`parameterSourceInfo$estimationInitialValues` must distinguish an asymmetric
policy. Collapsing it to "both roles use simulation starts" makes a run where
only the reference does read identically to one where neither does — the run
record then fails to state the choice actually made.

- [ ] **Step 7: Warn when a role policy cannot take effect**

`"simulation"` starts only differ from stored starts in `"rawres"` mode, but
silence is the wrong response to an explicitly-set option. Every other
mode-gated argument in `runSSEControl()` tracks `missing()` so an
explicitly-set-but-inapplicable value is never swallowed; these must too.

The two non-rawres modes are not alike, and the message must not pretend they
are. Under `"fixed"` a single generating vector is resolved once and reused, so
simulation starts are genuinely identical to model starts. Under `"covariance"`
each replicate draws a different vector, so simulation starts *would* differ —
the behaviour is unimplemented, not inert, and saying otherwise is misleading.

```r
  if (identical(referenceInitials, "simulation") && parameterSource != "rawres") {
    cli::cli_warn(c(
      "!" = "{.arg referenceInitials = \"simulation\"} has no effect unless {.arg parameterSource = \"rawres\"}.",
      "i" = "Under {.val fixed} every replicate shares one generating vector, so the two policies coincide.",
      "i" = "Under {.val covariance} the vectors differ per replicate, but simulation starts are not yet wired up."
    ))
  }
```

Do the same for `alternativeInitials`. A hard error is ruled out: it would break
the legitimate case of setting one role while the other stays at its default.

- [ ] **Step 7: Run the full suite and commit**

```bash
git add R/sse-control.R R/sse-helpers.R R/run-sse.R R/recompute-sse.R DESCRIPTION tests/testthat/test-control.R tests/testthat/test-run-sse.R
git commit -m "feat: role-specific estimation starting values"
```

---

## Task 5: Full-distribution PPE estimation

**Files:**
- Create: `R/sse-ppe.R`, `tests/testthat/test-sse-ppe.R`
- Modify: `R/plot-sse.R`

Uses `ppe_dofv()` and `fake_ppe_sse_object()` from
`tests/testthat/helper-ppe-fixtures.R`, written in Task 1 Step 2.

- [ ] **Step 1: Write the failing estimator tests**

```r
test_that("the MLE recovers a known noncentrality parameter", {
  for (spec in list(list(df = 1, ncp = 8), list(df = 3, ncp = 15), list(df = 5, ncp = 25))) {
    d <- ppe_dofv(n = 4000L, df = spec$df, ncp = spec$ncp, seed = 7L)
    fit <- .ppeChiSquareMle(d, df = spec$df)
    # 4000 draws gives an SE well under 5% of ncp for these settings
    expect_equal(fit$estimate, spec$ncp, tolerance = 0.05)
  }
})

test_that("the MLE recovers a known df when ncp is held at zero", {
  d <- ppe_dofv(n = 4000L, df = 4, ncp = 0, seed = 11L)
  fit <- .ppeChiSquareMle(d, ncp = 0)
  expect_equal(fit$estimate, 4, tolerance = 0.05)
})

test_that("the MLE agrees with a direct likelihood grid", {
  d <- ppe_dofv(n = 300L, df = 2, ncp = 6, seed = 13L)
  grid <- seq(0.1, 20, by = 0.01)
  ll <- vapply(grid, function(ncp) sum(stats::dchisq(d[d > 0], 2, ncp, log = TRUE)), numeric(1))

  expect_equal(.ppeChiSquareMle(d, df = 2)$estimate, grid[which.max(ll)], tolerance = 0.02)
})

test_that("non-positive values are counted, never hidden", {
  d <- c(-2, -1, 0, 3, 4, 5, 6)
  fit <- .ppeChiSquareMle(d, df = 1)

  expect_equal(fit$nRetained, 4L)
  expect_equal(fit$nNonPositive, 3L)
  expect_equal(fit$n, 7L)
})

test_that("a boundary solution is flagged rather than reported as a clean fit", {
  # retained mean below df drives the constrained MLE onto the lower bound
  d <- c(0.1, 0.2, 0.15, 0.3)
  fit <- .ppeChiSquareMle(d, df = 4)

  expect_true(fit$boundary)
  expect_lt(fit$estimate, 1e-6)
})

test_that("too few positive values is an error, not a silent fit", {
  expect_error(.ppeChiSquareMle(c(-1, 2), df = 1), "at least 2")
})
```

- [ ] **Step 2: Run to verify failure**

Expected: FAIL, `could not find function ".ppeChiSquareMle"`.

- [ ] **Step 3: Implement the estimator**

```r
.ppeChiSquareMle <- function(dofv, df = NULL, ncp = NULL) {
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
  estimate_ncp <- is.null(ncp)
  if (estimate_ncp == is.null(df)) {
    .abortSSE("Supply exactly one of {.arg df} and {.arg ncp}; the other is estimated.")
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
  fit <- stats::optim(par = init, fn = objective, lower = 1e-16,
                      method = "L-BFGS-B")
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
```

`init` may fall below `lower`; `optim()` projects it into the feasible region
and converges normally, so no guard is added. This was verified explicitly.

- [ ] **Step 4: Run to verify the estimator tests pass**

Expected: PASS.

- [ ] **Step 5: Write the failing `nonpositive` policy tests**

```r
test_that("the nonpositive policy controls warning and abort, never the counts", {
  d <- c(-1, -2, 3, 4, 5)

  expect_warning(.ppeApplyNonpositivePolicy(2L, 5L, "warn"), "2 of 5")
  expect_error(.ppeApplyNonpositivePolicy(2L, 5L, "error"), "2 of 5")
  expect_silent(.ppeApplyNonpositivePolicy(2L, 5L, "drop"))
})
```

- [ ] **Step 6: Implement the policy**

```r
.ppeApplyNonpositivePolicy <- function(nNonPositive, n, policy) {
  if (nNonPositive == 0L) return(invisible(NULL))
  msg <- c(
    "!" = "{nNonPositive} of {n} test statistic{?s} {?is/are} not positive.",
    "i" = "The noncentral chi-square has no support there, so they are excluded from the fit.",
    "i" = "The estimate is biased upward when this fraction is large."
  )
  switch(policy,
    error = .abortSSE(msg),
    warn = cli::cli_warn(msg),
    drop = invisible(NULL)
  )
}
```

- [ ] **Step 7: Add the method branch to `.ppePowerPlotData()`**

Resolve comparisons for the `distribution_mle` branch with `ppe = TRUE`, not
`FALSE`. The warning `.legacyComparisons()` emits when df is inferred exists
*specifically* because the MLE treats df as known — an inferred value silently
feeding a headline power number is precisely the failure this plan's decision
5 and 12 exist to prevent ("parameter-count differences remain a convenience
fallback, not an asserted truth"; "record every statistical choice"). That the
legacy `exceedance` path (via `.modelDegreesFreedom()` directly, not through
`.resolveComparisons()`) has always been silent here is not a reason for the
new default method to inherit the same gap — `exceedance` being silent is
exactly the legacy posture this plan replaces, not a precedent to match.

Also add a `df_source` column to the MLE output, copied from each resolved
comparison, so provenance survives even for a caller who suppresses the
warning. Verified empirically that its absence leaves literally no signal
anywhere — not a warning, not a column — that df was guessed rather than
declared.

Tests built on fixtures with no explicit `sseComparison` (which is most of this
package's existing test suite) will now see this warning under
`method = "distribution_mle"`. Wrap those specific calls in `expect_warning()`
or `suppressWarnings()` as appropriate — do not pass `ppe = FALSE` to silence
it package-wide.

Rename the legacy outputs so the two methods cannot be confused, per rigour
decision 9:

```r
  if (identical(method, "exceedance")) {
    rows[[row_index]]$threshold_exceedance_probability <- point_prob
    rows[[row_index]]$threshold_implied_ncp <- point_ncp
  }
```

- [ ] **Step 8: Verify the frozen contract still holds for `exceedance`**

Run: `NOT_CRAN=true Rscript -e 'pkgload::load_all("."); testthat::test_file("tests/testthat/test-sse-statistical-contract.R")'`
Expected: PASS with `method = "exceedance"` passed explicitly. Update the
contract test to name the method, since the default has changed.

- [ ] **Step 9: Commit**

```bash
git add R/sse-ppe.R tests/testthat/test-sse-ppe.R tests/testthat/helper-ppe-fixtures.R R/plot-sse.R tests/testthat/test-sse-statistical-contract.R
git commit -m "feat: distribution-based PPE noncentrality estimation"
```

---

## Task 6: Model-based PPE uncertainty

**Files:**
- Modify: `R/sse-ppe.R`, `tests/testthat/test-sse-ppe.R`, `R/plot-sse.R`

- [ ] **Step 1: Write the failing bootstrap tests**

```r
test_that("the bootstrap is reproducible under a fixed seed", {
  a <- .ppeParametricBootstrap(8, df = 1, nRetained = 100L,
                               bootstrapSamples = 50L, seed = 5L)
  b <- .ppeParametricBootstrap(8, df = 1, nRetained = 100L,
                               bootstrapSamples = 50L, seed = 5L)

  expect_equal(a$ci_lower, b$ci_lower)
  expect_equal(a$ci_upper, b$ci_upper)
})

test_that("the bootstrap interval brackets the point estimate", {
  res <- .ppeParametricBootstrap(8, df = 1, nRetained = 200L,
                                 bootstrapSamples = 200L, seed = 3L)

  expect_lt(res$ci_lower, 8)
  expect_gt(res$ci_upper, 8)
  expect_equal(res$n_successful, 200L)
  expect_equal(res$interval_type, "model_based")
})

test_that("the caller's RNG stream is left untouched", {
  set.seed(99L)
  before <- .Random.seed
  .ppeParametricBootstrap(8, df = 1, nRetained = 50L,
                          bootstrapSamples = 20L, seed = 5L)

  expect_equal(.Random.seed, before)
})

test_that("bootstrapSamples = 0 skips the bootstrap without changing the estimate", {
  res <- .ppeParametricBootstrap(8, df = 1, nRetained = 50L,
                                 bootstrapSamples = 0L, seed = 5L)

  expect_true(is.na(res$ci_lower))
  expect_equal(res$n_successful, 0L)
})

test_that("bootstrap failures are counted rather than aborting the fit", {
  res <- .ppeParametricBootstrap(1e-16, df = 1, nRetained = 3L,
                                 bootstrapSamples = 20L, seed = 5L)

  expect_equal(res$n_successful + res$n_failed, 20L)
})
```

- [ ] **Step 2: Run to verify failure**

Run: `NOT_CRAN=true Rscript -e 'pkgload::load_all("."); testthat::test_file("tests/testthat/test-sse-ppe.R")'`
Expected: FAIL, `could not find function ".ppeParametricBootstrap"`.

- [ ] **Step 3: Implement the seed guard and bootstrap**

```r
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

.ppeParametricBootstrap <- function(estimate, df, nRetained,
                                    bootstrapSamples = 1000L,
                                    target = c("ncp", "df"),
                                    conf.level = 0.95, seed = NULL) {
  target <- match.arg(target)
  empty <- list(ci_lower = NA_real_, ci_upper = NA_real_, n_successful = 0L,
                n_failed = 0L, interval_type = "model_based", draws = numeric(0))
  if (bootstrapSamples <= 0L) {
    return(empty)
  }

  draws <- .withPpeSeed(seed, {
    vapply(seq_len(bootstrapSamples), function(i) {
      tryCatch({
        if (identical(target, "ncp")) {
          .ppeChiSquareMle(stats::rchisq(nRetained, df = df, ncp = estimate),
                           df = df)$estimate
        } else {
          .ppeChiSquareMle(stats::rchisq(nRetained, df = estimate, ncp = 0),
                           ncp = 0)$estimate
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
    # Named model_based, never "empirical": rchisq() draws are strictly
    # positive, so no replicate reproduces the truncation the real data
    # underwent. This interval covers estimator variability under the fitted
    # model only, and never misspecification.
    interval_type = "model_based",
    draws = ok
  )
}
```

- [ ] **Step 4: Run to verify the bootstrap tests pass**

Run: `NOT_CRAN=true Rscript -e 'pkgload::load_all("."); testthat::test_file("tests/testthat/test-sse-ppe.R")'`
Expected: PASS.

- [ ] **Step 5: Write the failing ppeSummary test**

```r
test_that("ppeSummary reports estimate, interval, counts, and provenance", {
  sse <- fake_ppe_sse_object(df = 1, ncp = 10, n = 80L, seed = 21L)
  cmp <- sseComparison("simulation", "alt1", df = 1)

  s <- ppeSummary(sse, comparisons = cmp, bootstrapSamples = 100L, bootSeed = 4L)

  expect_equal(s$parameter, "ncp")
  expect_equal(s$mode, "power")
  expect_equal(s$df_source, "explicit")
  expect_equal(s$interval_type, "model_based")
  expect_true(all(c("n", "n_nonpositive", "n_bootstrap_successful",
                    "n_bootstrap_failed", "boundary", "probability") %in% names(s)))
  expect_gt(s$estimate, 0)
})

test_that("a comparison with a custom critical value is refused by PPE", {
  sse <- fake_ppe_sse_object(df = 1, ncp = 10, n = 40L, seed = 22L)
  cmp <- sseComparison("simulation", "alt1", criticalValue = 5)

  expect_error(ppeSummary(sse, comparisons = cmp), "noncentral chi-square")
})
```

- [ ] **Step 6: Implement .ppeFit() and ppeSummary()**

```r
.ppeDefaultSeed <- function(x, comparison) {
  # Derived from the run's own seed so a given run plots identically every
  # time while remaining distinct across runs.
  base <- as.integer(x$runInfo$seed %||% 1L)
  offset <- sum(utf8ToInt(comparison$label)) %% 100000L
  as.integer((base + offset) %% 2147483647L)
}

.ppeFit <- function(x, comparison, conf.level = 0.95, nonpositive = "warn",
                    bootstrapSamples = 1000L, bootSeed = NULL) {
  if (!isTRUE(comparison$ppeEligible)) {
    .abortSSE(c(
      "Comparison {.val {comparison$label}} has an explicit {.arg criticalValue} and no {.arg df}.",
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
```

- [ ] **Step 7: Render power curves and Type-I point-ranges**

A power comparison scales the noncentrality as `lambda_N = lambda * N / N_base`
and evaluates `pchisq(criticalValue, df, lambda_N, lower.tail = FALSE)`, with the
ribbon taken from the bootstrap limits rather than a binomial interval.

A Type-I comparison has no sample-size curve — the estimated quantity is `df`,
not a noncentrality that scales — so it renders instead as a point-range of the
estimated Type-I rate with its bootstrap interval, against a dashed nominal
`alpha` reference line.

```r
test_that("a type1 comparison renders a point-range, not a curve", {
  sse <- fake_ppe_sse_object(df = 1, ncp = 0, n = 200L, seed = 71L)
  p <- plotSSEPpePower(
    sse, comparisons = sseComparison("alt1", "simulation", df = 1),
    bootstrapSamples = 50L, bootSeed = 3L
  )

  geoms <- vapply(p$layers, function(l) class(l$geom)[[1L]], character(1))
  expect_true("GeomPointrange" %in% geoms)
  expect_false("GeomLine" %in% geoms)
})

test_that("mixing power and type1 comparisons warns rather than plotting both", {
  sse <- fake_ppe_sse_object(df = 1, ncp = 8, n = 100L, seed = 72L)
  expect_warning(
    plotSSEPpePower(sse, comparisons = list(
      sseComparison("simulation", "alt1", df = 1),
      sseComparison("alt1", "simulation", df = 1, label = "type1 check")
    ), bootstrapSamples = 20L, bootSeed = 3L),
    "omitted"
  )
})
```

Implement by branching on `unique(mode)` across the resolved comparisons. When
both appear, render the power panel and warn that the Type-I comparisons were
omitted, directing the reader to `ppeSummary()`.

- [ ] **Step 8: Run the full suite and commit**

Run: `NOT_CRAN=true Rscript -e 'pkgload::load_all("."); testthat::test_dir("tests/testthat", reporter="progress", stop_on_failure=FALSE)'`
Expected: `FAIL 0`.

```bash
git add R/sse-ppe.R tests/testthat/test-sse-ppe.R R/plot-sse.R NAMESPACE
git commit -m "feat: parametric bootstrap uncertainty, ppeSummary(), and type1 rendering"
```

---

## Task 7: PPE distribution diagnostics

**Files:**
- Create: `R/plot-sse-ppe-diagnostics.R`, `tests/testthat/test-plot-sse-ppe-diagnostics.R`
- Modify: `NAMESPACE`; create `_pkgdown.yml`

- [ ] **Step 1: Write the failing discrepancy tests**

```r
test_that("Cramer-von Mises separates a correct model from a misspecified one", {
  good <- ppe_dofv(n = 400L, df = 1, ncp = 9, seed = 31L)
  bad <- c(ppe_dofv(200L, df = 1, ncp = 2, seed = 32L),
           ppe_dofv(200L, df = 1, ncp = 30, seed = 33L))

  cvm_good <- .ppeCramerVonMises(good[good > 0], df = 1, ncp = 9)
  cvm_bad <- .ppeCramerVonMises(bad[bad > 0], df = 1,
                                ncp = .ppeChiSquareMle(bad, df = 1)$estimate)

  expect_lt(cvm_good, cvm_bad)
})

test_that("the diagnostic p-value is small for a misspecified mixture", {
  bad <- c(ppe_dofv(200L, df = 1, ncp = 2, seed = 34L),
           ppe_dofv(200L, df = 1, ncp = 30, seed = 35L))
  fit <- .ppeChiSquareMle(bad, df = 1)

  p <- .ppeDiagnosticPValue(bad[bad > 0], df = 1, estimate = fit$estimate,
                            bootstrapSamples = 200L, seed = 36L)

  expect_lt(p, 0.05)
})

test_that("the diagnostic p-value is not small when the model is correct", {
  good <- ppe_dofv(n = 400L, df = 1, ncp = 9, seed = 37L)
  fit <- .ppeChiSquareMle(good, df = 1)

  p <- .ppeDiagnosticPValue(good[good > 0], df = 1, estimate = fit$estimate,
                            bootstrapSamples = 200L, seed = 38L)

  expect_gt(p, 0.05)
})
```

- [ ] **Step 2: Run to verify failure**

Expected: FAIL, `could not find function ".ppeCramerVonMises"`.

- [ ] **Step 3: Implement the discrepancy statistic and its p-value**

```r
.ppeCramerVonMises <- function(x, df, ncp) {
  x <- sort(x)
  n <- length(x)
  u <- stats::pchisq(x, df = df, ncp = ncp)
  1 / (12 * n) + sum((u - (2 * seq_len(n) - 1) / (2 * n))^2)
}

.ppeDiagnosticPValue <- function(x, df, estimate, bootstrapSamples = 1000L,
                                 seed = NULL) {
  observed <- .ppeCramerVonMises(x, df = df, ncp = estimate)
  n <- length(x)
  null_stats <- .withPpeSeed(seed, {
    vapply(seq_len(bootstrapSamples), function(i) {
      draw <- stats::rchisq(n, df = df, ncp = estimate)
      refit <- tryCatch(.ppeChiSquareMle(draw, df = df)$estimate,
                        error = function(e) NA_real_)
      if (is.na(refit)) {
        return(NA_real_)
      }
      .ppeCramerVonMises(draw[draw > 0], df = df, ncp = refit)
    }, numeric(1))
  })
  null_stats <- null_stats[is.finite(null_stats)]
  if (length(null_stats) == 0L) {
    return(NA_real_)
  }
  mean(null_stats >= observed)
}
```

- [ ] **Step 4: Run to verify the discrepancy tests pass**

Expected: PASS.

- [ ] **Step 5: Write the failing plot tests**

```r
test_that("plotSSEPpeDiagnostics returns a ggplot carrying auditable data", {
  sse <- fake_ppe_sse_object(df = 1, ncp = 10, n = 80L, seed = 41L)
  p <- plotSSEPpeDiagnostics(
    sse, comparisons = sseComparison("simulation", "alt1", df = 1),
    bootstrapSamples = 50L, bootSeed = 2L
  )

  expect_s3_class(p, "ggplot")
  diag <- attr(p, "ppeDiagnostics")
  expect_true(all(c("comparison", "cvm", "p_value", "n", "n_nonpositive") %in%
                    names(diag)))
})

test_that("the diagnostic subtitle names the excluded count", {
  sse <- fake_ppe_sse_object(df = 1, ncp = 10, n = 80L, seed = 42L,
                             nNonPositive = 6L)
  p <- suppressWarnings(plotSSEPpeDiagnostics(
    sse, comparisons = sseComparison("simulation", "alt1", df = 1),
    bootstrapSamples = 20L, bootSeed = 2L
  ))

  expect_match(p$labels$subtitle, "6")
})
```

- [ ] **Step 6: Implement plotSSEPpeDiagnostics()**

Build the ECDF panel by overlaying `ggplot2::stat_ecdf()` of the retained
statistics on the fitted `stats::pchisq()` curve, with a pointwise envelope
formed from the bootstrap draws. Add a quantile-quantile panel against an
identity line. Set the subtitle from `n`, `n_nonpositive`, `df`, and
`df_source`, and the caption to state that the envelope is pointwise rather
than simultaneous. Attach the diagnostic table:

```r
  attr(p, "ppeDiagnostics") <- data.frame(
    comparison = cmp$label, cvm = cvm, p_value = p_value,
    n = fit$n, n_nonpositive = fit$nNonPositive,
    df = cmp$df, df_source = cmp$dfSource,
    stringsAsFactors = FALSE
  )
```

Diagnostics are evidence about approximation adequacy, not a pass/fail
certification; the documentation must say so.

- [ ] **Step 7: Run and commit**

```bash
git add R/plot-sse-ppe-diagnostics.R tests/testthat/test-plot-sse-ppe-diagnostics.R NAMESPACE _pkgdown.yml
git commit -m "feat: PPE distribution adequacy diagnostics"
```

---

## Task 8: Validate sample-size proportionality

**Files:**
- Create: `R/sse-ppe-validation.R`, `tests/testthat/test-sse-ppe-validation.R`
- Modify: `NAMESPACE`, `_pkgdown.yml`

- [ ] **Step 1: Write the failing tests**

```r
test_that("proportional fixtures give a consistent lambda per subject", {
  runs <- list(
    fake_ppe_sse_object(df = 1, ncp = 5,  n = 400L, subjects = 20L, seed = 51L),
    fake_ppe_sse_object(df = 1, ncp = 10, n = 400L, subjects = 40L, seed = 52L),
    fake_ppe_sse_object(df = 1, ncp = 20, n = 400L, subjects = 80L, seed = 53L)
  )
  v <- validateSSEPpeScaling(
    runs, comparisons = sseComparison("simulation", "alt1", df = 1)
  )

  expect_equal(v$table$lambda_per_subject, rep(0.25, 3), tolerance = 0.15)
  expect_false(v$nonlinear)
})

test_that("a deliberately nonlinear fixture is flagged", {
  runs <- list(
    fake_ppe_sse_object(df = 1, ncp = 5, n = 400L, subjects = 20L, seed = 54L),
    fake_ppe_sse_object(df = 1, ncp = 8, n = 400L, subjects = 40L, seed = 55L),
    fake_ppe_sse_object(df = 1, ncp = 9, n = 400L, subjects = 80L, seed = 56L)
  )
  expect_warning(
    v <- validateSSEPpeScaling(
      runs, comparisons = sseComparison("simulation", "alt1", df = 1)
    ),
    "not proportional"
  )
  expect_true(v$nonlinear)
})

test_that("two study sizes give only a descriptive ratio", {
  runs <- list(
    fake_ppe_sse_object(df = 1, ncp = 5,  n = 200L, subjects = 20L, seed = 57L),
    fake_ppe_sse_object(df = 1, ncp = 10, n = 200L, subjects = 40L, seed = 58L)
  )
  v <- validateSSEPpeScaling(
    runs, comparisons = sseComparison("simulation", "alt1", df = 1)
  )

  expect_true(is.na(v$lackOfFit))
  expect_equal(nrow(v$table), 2L)
})

test_that("incompatible runs are refused rather than silently pooled", {
  runs <- list(
    fake_ppe_sse_object(df = 1, ncp = 5, n = 200L, subjects = 20L, seed = 59L),
    fake_ppe_sse_object(df = 2, ncp = 5, n = 200L, subjects = 40L, seed = 60L)
  )
  expect_error(
    validateSSEPpeScaling(
      runs, comparisons = sseComparison("simulation", "alt1", df = 1)
    ),
    "incompatible"
  )
})
```

- [ ] **Step 2: Run to verify failure**

Expected: FAIL, `could not find function "validateSSEPpeScaling"`.

- [ ] **Step 3: Implement the validator**

```r
#' Validate the proportional-noncentrality assumption across study sizes
#' @export
validateSSEPpeScaling <- function(runs, comparisons = NULL, conf.level = 0.95,
                                  bootstrapSamples = 1000L, tolerance = 0.25) {
  if (!is.list(runs) || length(runs) < 2L) {
    .abortSSE("{.arg runs} must be a list of at least 2 completed SSE objects.")
  }
  lapply(runs, .assertSSEObject)

  fits <- lapply(runs, function(run) {
    cmp <- .resolveComparisons(run, comparisons %||% run$runInfo$comparisons,
                               ppe = TRUE)[[1L]]
    list(cmp = cmp,
         fit = .ppeFit(run, cmp, conf.level = conf.level,
                       bootstrapSamples = bootstrapSamples),
         subjects = .studySampleSize(run)$size)
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

  lack_of_fit <- NA_real_
  if (nrow(table) >= 3L) {
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
```

- [ ] **Step 4: Run and commit**

```bash
git add R/sse-ppe-validation.R tests/testthat/test-sse-ppe-validation.R NAMESPACE _pkgdown.yml
git commit -m "feat: validate the proportional noncentrality assumption"
```

---

## Task 9: Correct Monte Carlo uncertainty for parameter summaries

**Files:**
- Modify: `R/sse-helpers.R`, `R/sse-methods.R`
- Create: `tests/testthat/test-sse-parameter-summary.R`

- [ ] **Step 1: Write the failing MCSE tests**

```r
test_that("bias MCSE is the SD of replicate errors over sqrt(n)", {
  errors <- c(0.1, -0.2, 0.3, 0.05)
  expect_equal(.mcseFromErrors(errors), stats::sd(errors) / sqrt(4))
})

test_that("MCSE is nonnegative for a negative fixed truth", {
  s <- .parameterSummaryRow(estimates = c(-1.8, -2.2, -2.1), truth = -2)

  expect_gte(s$mcse_bias, 0)
  expect_gte(s$mcse_relative_bias, 0)
})

test_that("zero truths are excluded from relative bias", {
  s <- .parameterSummaryRow(estimates = c(0.1, 0.2), truth = 0)

  expect_true(is.na(s$mcse_relative_bias))
  expect_equal(s$n_effective_relative, 0L)
})

test_that("a single replicate yields NA intervals with its count retained", {
  s <- .parameterSummaryRow(estimates = 1.5, truth = 1)

  expect_true(is.na(s$ci_bias_lower))
  expect_equal(s$n_effective, 1L)
})
```

- [ ] **Step 2: Run to verify failure**

Expected: FAIL, `could not find function ".mcseFromErrors"`.

- [ ] **Step 3: Implement the replicate-level statistics**

```r
.mcseFromErrors <- function(errors) {
  errors <- errors[is.finite(errors)]
  n <- length(errors)
  if (n < 2L) {
    return(NA_real_)
  }
  # An SD is nonnegative by construction, so this is nonnegative even when the
  # truth is negative -- which the previous `rse` field was not.
  stats::sd(errors) / sqrt(n)
}
```

`.parameterSummaryRow()` builds absolute errors as `estimate - truth` and
relative errors as `(estimate - truth) / truth`, dropping non-finite values and
excluding zero truths from the relative branch. Normal-approximation intervals
are emitted only when `n_effective >= 2` and the MCSE is finite and positive;
otherwise the interval bounds are `NA_real_` and the count is still reported.

- [ ] **Step 4: Add the new tidy-summary fields**

Add `mcse_bias`, `ci_bias_lower`, `ci_bias_upper`, `mcse_relative_bias`,
`ci_relative_bias_lower`, `ci_relative_bias_upper`, `n_effective`, and
`n_effective_relative`.

Replace the ambiguous `rse` with `mcse_relative_bias`, retaining `rse` as a
documented superseded alias for one release. It is a table field rather than a
callable argument, so `lifecycle` cannot warn when it is read; document that
plainly rather than implying a warning exists.

- [ ] **Step 5: Run and commit**

```bash
git add R/sse-helpers.R R/sse-methods.R tests/testthat/test-sse-parameter-summary.R tests/testthat/test-sse-statistical-contract.R
git commit -m "feat: correctly named Monte Carlo uncertainty for parameter summaries"
```

---

## Task 10: Parameter-draw adequacy summaries

**Files:**
- Create: `R/sse-parameter-diagnostics.R`, `tests/testthat/test-sse-parameter-diagnostics.R`
- Modify: `R/run-sse.R`, `R/sse-methods.R`, `NAMESPACE`, `tests/testthat/helper-ppe-fixtures.R`, `_pkgdown.yml`

- [ ] **Step 1: Add the parameter-draw fixture**

`parameterDrawSummary()` reads the realised generating values from
`sse$initialValues` and the intended targets from
`runInfo$parameterSourceInfo`. This fixture synthesises both, so the summary can
be tested without running a real covariance-mode SSE.

```r
# tests/testthat/helper-ppe-fixtures.R

# Draws OMEGA the way each covariance mode does, so the known artefacts are
# reproduced exactly: log-Cholesky inflates the raw-scale mean by
# exp(se^2 / (2 * omega^2)); the inverse-Wishart route does not.
fake_draw_sse_object <- function(mode = c("joint", "independent_iw"),
                                 omega = 0.6, se = 0.3, n = 200L,
                                 thetaLower = -Inf, seed = 61L) {
  mode <- match.arg(mode)
  has <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old <- if (has) get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit({
    if (has) assign(".Random.seed", old, envir = globalenv())
  }, add = TRUE)
  set.seed(seed)

  omega_draws <- if (identical(mode, "joint")) {
    exp(stats::rnorm(n, mean = log(omega), sd = se / omega))
  } else {
    nu <- 3 + 2 * (omega / se)^2
    vapply(seq_len(n), function(i) {
      drop(rxode2::cvPost(nu, matrix(omega * (nu - 2) / nu, 1L, 1L), n = 1L))
    }, numeric(1))
  }
  theta_draws <- stats::rnorm(n, mean = 0.45, sd = 0.2)

  sse <- fake_ppe_sse_object(n = n, seed = seed)
  sse$initialValues <- data.frame(
    replicate = rep(seq_len(n), 2L),
    parameter = rep(c("tka", "omega(eta.ka,eta.ka)"), each = n),
    value = c(theta_draws, omega_draws),
    stringsAsFactors = FALSE
  )
  sse$runInfo$parameterSource <- "covariance"
  sse$runInfo$covarianceDraw <- mode
  sse$runInfo$parameterSourceInfo <- list(
    covarianceDraw = mode,
    targets = data.frame(
      parameter = c("tka", "omega(eta.ka,eta.ka)"),
      target_mean = c(0.45, omega),
      target_sd = c(0.2, se),
      lower = c(thetaLower, 0),
      stringsAsFactors = FALSE
    )
  )
  sse
}
```

- [ ] **Step 2: Write the failing tests**

```r
test_that("joint mode exposes raw-scale OMEGA mean drift", {
  sse <- fake_draw_sse_object(mode = "joint", omega = 0.6, se = 0.3,
                              n = 500L, seed = 61L)
  s <- parameterDrawSummary(sse)
  row <- s[s$parameter == "omega(eta.ka,eta.ka)", ]

  # The log-Cholesky back-transformation inflates the raw mean by
  # exp(SE^2 / (2 * Omega^2)); at SE/Omega = 0.5 that is about 13%.
  expect_gt(row$realized_mean / row$target_mean, 1.05)
  expect_true(row$mean_drift_flag)
})

test_that("independent_iw reports binding degrees of freedom per block", {
  sse <- fake_draw_sse_object(mode = "independent_iw", omega = 0.6, se = 0.3,
                              n = 200L, seed = 62L)
  s <- parameterDrawSummary(sse)

  expect_true(all(s$binding_nu[!is.na(s$binding_nu)] > 0))
})

test_that("every drawn OMEGA block is positive definite", {
  sse <- fake_draw_sse_object(mode = "joint", omega = 0.6, se = 0.3,
                              n = 200L, seed = 63L)
  s <- parameterDrawSummary(sse)

  expect_equal(sum(s$n_not_positive_definite, na.rm = TRUE), 0L)
})

test_that("out-of-bound THETA draws are counted, never resampled", {
  sse <- fake_draw_sse_object(mode = "joint", thetaLower = 0.5,
                              n = 200L, seed = 64L)
  s <- parameterDrawSummary(sse)

  expect_gt(sum(s$n_out_of_domain, na.rm = TRUE), 0L)
  expect_equal(nrow(s), length(unique(s$parameter)))
})
```

- [ ] **Step 3: Run to verify failure**

Expected: FAIL, `could not find function "parameterDrawSummary"`.

- [ ] **Step 4: Implement the summary**

Report realised mean, SD, median, central quantiles, minimum, maximum, and
finite count for every generating parameter. During `runSSE()`, persist enough
target metadata for covariance mode to compare realised THETA means and SDs and
OMEGA means with the intended local approximation.

For `joint`, report raw-scale OMEGA mean drift and empirical dependence for
covered parameters, without claiming exact covariance preservation. For
`independent_iw`, report each block's binding degrees of freedom and intended
versus realised diagonal dispersion. Check every generated OMEGA block for
positive definiteness. Compare THETA draws against recoverable finite model
bounds and count out-of-domain draws, never truncating or resampling.

Store the table in the returned object and in `sse_summary.rds`. Warn only on
prespecified actionable conditions; always retain the full table.

- [ ] **Step 5: Run and commit**

```bash
git add R/sse-parameter-diagnostics.R tests/testthat/test-sse-parameter-diagnostics.R R/run-sse.R R/sse-methods.R NAMESPACE _pkgdown.yml
git commit -m "feat: parameter-draw adequacy summaries"
```

---

## Task 11: Documentation, migration, and package validation

**Files:**
- Modify: `docs/sse-technical-reference.md`, `README.md`, `vignettes/runSSE.Rmd`, `NEWS.md`, `_pkgdown.yml`, roxygen comments

The technical reference is updated **last**, in Step 5, once everything else is
implemented and verified. It is a description of what the package does, so
writing it before the validation sweep passes risks documenting intent rather
than behaviour.

- [ ] **Step 1: Write the NEWS entries**

```markdown
## Breaking changes

* `plotSSEPpePower()` now defaults to `method = "distribution_mle"`, which
  estimates one noncentrality parameter per comparison by maximum likelihood
  from the whole test-statistic distribution. Previous releases inverted the
  empirical exceedance probability separately at each threshold, which could
  imply a different effect size at every threshold. Plots and derived numbers
  will change. Pass `method = "exceedance"` to restore the previous estimator;
  its outputs are now named `threshold_exceedance_probability` and
  `threshold_implied_ncp` so the two cannot be confused.
* Parameter summaries replace `rse` with `mcse_relative_bias`, which is always
  nonnegative and defined for varying truths. `rse` remains as a superseded
  alias for one release.
```

- [ ] **Step 2: Document the limitations honestly**

In `?plotSSEPpePower`, under a **Limitations** heading, state that the fit
excludes non-positive statistics the noncentral chi-square cannot produce and is
therefore biased upward when they are common; that the parametric bootstrap
covers estimator variability under the fitted model and never model
misspecification; and that proportional scaling is an assumption testable only
with `validateSSEPpeScaling()`.

State that an explicit comparison does not prove nesting or justify its null
distribution, and that boundary hypotheses need user-supplied critical values or
empirical null calibration.

- [ ] **Step 3: Regenerate documentation with the pinned toolchain**

Run: `Rscript -e 'devtools::document()'`
Expected: `Config/roxygen2/version: 8.1.0` unchanged. Compare the `NAMESPACE`
diff as a sorted set and confirm only the intended exports were added:
`sseComparison`, `comparisonSummary`, `ppeSummary`, `plotSSEPpeDiagnostics`,
`validateSSEPpeScaling`, `parameterDrawSummary`, `plotSSEParameterDraws`.

- [ ] **Step 4: Run the full validation sweep**

```bash
NOT_CRAN=true Rscript -e 'pkgload::load_all("."); testthat::test_dir("tests/testthat", reporter="progress", stop_on_failure=FALSE)'
Rscript -e 'devtools::check()'
```

Expected: `FAIL 0`; `R CMD check` reports 0 errors and 0 warnings.

**`R CMD check` is not clean today.** Measured on a clean `git archive` of the
tree at Task 3, so these are pre-existing and must be cleared here rather than
discovered at the end:

| Finding | Severity | Fix |
| --- | --- | --- |
| `plot-nlmixr2SSE.Rd` documents an argument `object` that is absent from `\usage` | **WARNING** | Reconcile the roxygen block against the actual formals, then re-document. This alone fails a `error_on = "warning"` gate. |
| `Namespace in Imports field not imported from: 'tidyr'` | NOTE | Either use it or drop it from `DESCRIPTION`. |
| `no visible binding for global variable` for `delta_ofv`, `model_label`, `value`, `role`, `estimate` | NOTE | Standard ggplot2 non-standard-evaluation artefact. Declare them via `utils::globalVariables()` in a package-level file. Tasks 7 and 10 add plots and will add more names — collect them all here. |

Two further vignette WARNINGs (`Files in the 'vignettes' directory but no files
in 'inst/doc'`) appear only when checking with `--no-build-vignettes`; they are
an artefact of that flag, not a real defect. Run the final check WITHOUT it.

Note also, because it is a natural false alarm: `stats::` and `utils::` are used
throughout `R/` without appearing in `DESCRIPTION`, and that is fine. `R CMD
check` does not require base-priority packages to be declared, and it does not
flag them. Do not "fix" this.

For test dependencies the rule is different: a bare `::` into a package in
neither `Imports` nor `Suggests` produces
`checking for unstated dependencies in 'tests' ... WARNING`. Either declare the
package or use base R instead.

- [ ] **Step 5: Update the technical reference**

This is the final substantive step. Do it only once Step 4 is green, so the
guide describes what the package does rather than what the plan intended.

`docs/sse-technical-reference.md` may carry uncommitted edits. Run
`git diff docs/sse-technical-reference.md` first and build on them; do not
overwrite. Section line numbers below are from the 734-line version at the time
of planning and will have moved.

**`## OFV differences, empirical power, and Type I error`** — the internal
difference is currently defined as
\(\Delta_i=\operatorname{OFV}_{\mathrm{reference},i}-\operatorname{OFV}_{\mathrm{alternative},i}\)
with `direction` labels encoding tail directions only. Replace this with the
comparison-based definition: the statistic is
\(T_i=\operatorname{OFV}_{\mathrm{reduced},i}-\operatorname{OFV}_{\mathrm{full},i}\),
computed from an explicit `sseComparison()` rather than inferred from a sign.
State that the mode (`power` or `type1`) follows from which member was the
simulation model, and that a comparison naming neither is refused because no
hypothesis is then known true. Document the paired-evaluable denominator and
every reported count: `n_attempted`, `n_full_evaluable`, `n_reduced_evaluable`,
`n_paired_evaluable`, `n_excluded`, `n_exceeding`.

**`## Parametric power estimation`** — rewrite with the distribution-based
method as primary. It must state:

- the estimand: one \(\lambda\) per comparison, fitted by maximum likelihood to
  the whole distribution of positive \(T_i\), replacing the threshold-specific
  inversion described in the existing paragraph beginning "The inversion is
  performed independently for every threshold";
- that the fit retains only \(T_i>0\) while using the *unconditional* noncentral
  chi-square density, so it is biased upward, and that the excluded count is
  always reported;
- that the constrained MLE sits at the lower bound when the retained mean falls
  below `df`, giving power exactly equal to `alpha` and a degenerate interval —
  measured at **31/200** replicates for ncp 0.5, df 4, n 60, so this is common
  rather than exotic;
- that the parametric bootstrap interval is `model_based`: it draws from
  `rchisq()`, which is strictly positive, so no replicate reproduces the
  truncation the real data underwent, and it never covers misspecification;
- that \(\lambda(n)=\lambda_0 n/n_0\) is an assumption, testable only by
  `validateSSEPpeScaling()` across directly simulated study sizes;
- what the ECDF, QQ, and Cramér–von Mises diagnostics can and cannot show —
  evidence about approximation adequacy, never a pass/fail certification;
- that `method = "exceedance"` retains the previous estimator verbatim, with
  outputs renamed `threshold_exceedance_probability` and
  `threshold_implied_ncp`.

Delete the sentence stating that `plotSSEPpePower()` is unavailable when
`estimateSimulation = FALSE`; Task 5 removes that restriction.

Replace the degrees-of-freedom paragraph. \(d=\max\{1,p_{sim}-p_{alt}\}\) is now
the fallback, not the definition: df comes from `sseComparison(df=)`, and the
parameter-count route is marked `dfSource = "parameter_count"` and warns when
PPE consumes it.

**`## Parameter summaries`** — replace `rse` with `mcse_relative_bias` and
document `mcse_bias`, both `ci_*` pairs, `n_effective`, and
`n_effective_relative`. State that MCSE is computed from replicate-level errors
and is therefore nonnegative even for a negative fixed truth, which `rse` was
not.

**`## Public interface`** — add `sseComparison()`, `comparisonSummary()`,
`ppeSummary()`, `plotSSEPpeDiagnostics()`, `validateSSEPpeScaling()`,
`parameterDrawSummary()`, and `plotSSEParameterDraws()`.

**`## Statistically material differences from PsN`** — this section already
carries hand-written edits. Amend rather than rewrite: the PPE estimator is no
longer a material difference, since both now fit the Ueckert (2016)
distribution-based noncentrality. What remains different is that `nlmixr2sse`
derives the comparison mode structurally from which model was simulated, where
PsN reads it from a model filename suffix (`_base|_r|_red|_reduced`); and that
`nlmixr2sse` ships the ECDF diagnostic PsN's templates also provide, so the
existing sentence saying it "does not implement that diagnostic" must go.

**`## Interpretation checklist and limitations`** — the existing bullets on
multiplicity and starting values already anticipate this work. Add bullets for
the paired-evaluable denominator, the non-positive \(T_i\) count, and whether
the proportional-scaling assumption was checked.

- [ ] **Step 6: Verify the guide against the code**

For each numbered claim added in Step 5, confirm the behaviour exists:

```bash
Rscript -e 'pkgload::load_all("."); print(names(formals(plotSSEPpePower)))'
Rscript -e 'pkgload::load_all("."); print(names(formals(sseComparison)))'
Rscript -e 'pkgload::load_all("."); print(names(ppeSummary(fake_ppe_sse_object())))'
```

Expected: every argument and column named in the guide appears. A guide claim
with no corresponding output is a defect in one or the other; resolve it rather
than softening the wording.

- [ ] **Step 7: Commit**

```bash
git add docs/sse-technical-reference.md README.md vignettes/runSSE.Rmd NEWS.md _pkgdown.yml man NAMESPACE
git commit -m "docs: document distribution-based PPE and Monte Carlo uncertainty"
```

---

## Statistical validation matrix

| Concern | Required validation | Task |
| --- | --- | --- |
| OFV direction | Hand-calculated full/reduced fixture with known sign | 1, 3 |
| Paired denominator | Asymmetric missingness and filtering fixture | 3 |
| Empirical power CI | Exact binomial calculation from reported counts | 3 |
| Initialization sensitivity | Same datasets and truths, different fit starts | 4 |
| PPE NCP MLE | Direct likelihood grid and simulated known-NCP fixtures | 5 |
| PPE df MLE | Known-df fixture with ncp held at zero | 5 |
| Boundary solutions | Retained mean below df drives the estimate to the bound | 5 |
| PPE bootstrap | Fixed-seed reproducibility, RNG restoration, failure accounting | 6 |
| PPE distribution diagnostic | Correct-model and misspecified-mixture fixtures | 7 |
| Sample-size scaling | Linear and deliberately nonlinear multi-size fixtures | 8 |
| Parameter bias MCSE | Hand calculation for fixed, negative, and varying truths | 9 |
| `independent_iw` draws | Mean, binding/nonbinding SD, positive-definiteness | 10 |
| `joint` draws | Positive-definiteness, mean drift, empirical dependence | 10 |
| Resume/recompute | Statistical definitions and RNG provenance unchanged | 2, 4 |

Monte Carlo tolerances must be chosen before looking at results, justified from
the expected Monte Carlo standard error, and wide enough for every supported
platform — never widened afterwards to make a test pass.

## Risks and mitigations

**Distribution-MLE PPE can look authoritative when its assumptions fail.**
Require explicit comparisons, report non-positive statistics, ship ECDF and QQ
diagnostics with a discrepancy p-value, and support direct multi-size
validation. Never describe a bootstrap band as covering misspecification.

**Changing the default estimator silently changes historical results.** The
`exceedance` method is retained verbatim under an explicit name, its output
fields are renamed so the two cannot be confused, method metadata is attached to
every output, and the change is announced in NEWS as breaking.

**Excluding failed fits can bias operating characteristics.** Always expose
attempted, evaluable, and paired counts; warn below a configurable paired
fraction; never impute.

**Starting at the generating truth can overstate numerical performance.** Keep
model starts as the default, record starts by role, and make same-dataset
sensitivity comparisons straightforward.

**Parameter-count degrees of freedom can be wrong.** Mark inferred df, warn when
PPE consumes it, and make explicit comparisons the documented workflow. Boundary
and constrained hypotheses require user justification or empirical null
calibration.

**Parametric bootstrap is costly.** Bootstrap only the fitted reference
distribution, use a keyed reproducible seed, and allow `bootstrapSamples = 0`.

## Out of scope

- NONMEM control-stream, `$PRIOR`, NWPRI, TNPRI, MSF, and covariance-file support.
- Exact PsN random-number or output-file reproduction, and PsN's PDF layout.
- Natural-scale multivariate-Normal OMEGA draws that can be indefinite.
- Automatic proof that models are nested or that their OFVs are comparable.
- Automatic selection of boundary-mixture null distributions.
- Treating covariance approximations as Bayesian posteriors.
- Silently clipping, truncating, repairing, or resampling invalid parameter draws.
- Changes to the underlying `nlmixr2est` fitting algorithms.
- Fitting a renormalised conditional likelihood to correct the truncation bias.
  That is a genuine improvement on the published method and needs its own spec
  and validation, not a silent divergence.

## Definition of done

1. Every reported comparison has an auditable definition and a paired denominator.
2. Empirical probabilities and parameter biases carry correctly named Monte
   Carlo uncertainty, distinct from model-based intervals.
3. Starting-value policy is role-specific, reproducible, and separable from the
   simulated-data distribution.
4. Distribution-MLE PPE estimates one noncentrality parameter per comparison,
   with bootstrap uncertainty and distributional diagnostics.
5. Multi-size runs can assess the proportional-noncentrality assumption.
6. Realized parameter draws have standard adequacy diagnostics.
7. Legacy outputs remain reproducible through explicit compatibility options.
8. Documentation states limitations without claiming the software can verify
   assumptions the user supplied.

## References

1. Ueckert S, Karlsson MO, Hooker AC. Accelerating Monte Carlo power studies
   through parametric power estimation. *J Pharmacokinet Pharmacodyn*.
   2016;43:223-234. [doi:10.1007/s10928-016-9468-y](https://doi.org/10.1007/s10928-016-9468-y).
2. Self SG, Liang K-Y. Asymptotic properties of maximum likelihood estimators
   and likelihood ratio tests under nonstandard conditions. *J Am Stat Assoc*.
   1987;82:605-610. [doi:10.1080/01621459.1987.10478472](https://doi.org/10.1080/01621459.1987.10478472).
3. Oehlert GW. A note on the delta method. *Am Stat*. 1992;46:27-29.
   [doi:10.1080/00031305.1992.10475842](https://doi.org/10.1080/00031305.1992.10475842).
4. Pinheiro JC, Bates DM. Unconstrained parametrizations for
   variance-covariance matrices. *Stat Comput*. 1996;6:289-296.
   [doi:10.1007/BF00140873](https://doi.org/10.1007/BF00140873).
5. Dosne A-G, Bergstrand M, Karlsson MO. An automated sampling importance
   resampling procedure for estimating parameter uncertainty. *J Pharmacokinet
   Pharmacodyn*. 2017;44:509-520. [doi:10.1007/s10928-017-9542-0](https://doi.org/10.1007/s10928-017-9542-0).
