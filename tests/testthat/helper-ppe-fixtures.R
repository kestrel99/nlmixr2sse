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
                                   subjects = 12L,
                                   seed = 42L) {
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
