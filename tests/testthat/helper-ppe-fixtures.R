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

# A Type-I-flavoured counterpart to fake_ppe_sse_object(). That fixture fixes
# the SIMULATION model's OFV low and lets the alternative float above it by
# the draw -- valid for POWER comparisons, where the simulation model is
# cmp$full (the richer, TRUE model that always fits at least as well). A
# Type-I comparison instead has cmp$reduced == simulation: the alternative is
# cmp$full, so a valid (non-negative) nested test statistic needs the
# ALTERNATIVE model's OFV to be the lower one, with the simulation model's
# floating above it by a CENTRAL (ncp = 0) chi-square draw -- the opposite
# assignment from fake_ppe_sse_object(). Reusing that fixture's sim-low
# convention for a Type-I comparison instead produces test statistics that
# are negative with probability 1, which .ppeChiSquareMle() correctly refuses
# to fit ("needs at least 2 positive test statistics") -- not a bug in the
# estimator, just the wrong fixture for the job.
fake_ppe_type1_sse_object <- function(df = 1, n = 200L, seed = 101L, subjects = 12L) {
  d <- ppe_dofv(n = n, df = df, ncp = 0, seed = seed)
  alt_floor <- rep(1000, n)
  fake_paired_sse_object(full_ofv = alt_floor + d, reduced_ofv = alt_floor,
                         subjects = subjects, seed = seed)
}

# A THREE-model fixture -- simulation, "alt1", "alt2" -- built so a single
# call can genuinely diagnose one power comparison (simulation vs alt1) and
# one Type-I comparison (alt2 vs simulation) together. Neither
# fake_ppe_sse_object() nor fake_ppe_type1_sse_object() can do this alone:
# each has only two models, and the OFV sign convention that makes one
# comparison direction valid (positive test statistics) makes the reversed
# direction on that SAME pair negative with probability 1 (see
# fake_ppe_type1_sse_object()'s header). Here the simulation model's OFV is
# fixed, alt1's floats ABOVE it by a noncentral chi-square draw (valid power
# direction: simulation vs alt1), and alt2's floats BELOW it by a central
# chi-square draw (valid Type-I direction: alt2 vs simulation) -- both hold
# simultaneously because they involve different model pairs.
fake_ppe_mixed_sse_object <- function(df = 1, powerNcp = 10, n = 80L, seed = 101L,
                                      subjects = 12L) {
  power_d <- ppe_dofv(n = n, df = df, ncp = powerNcp, seed = seed)
  type1_d <- ppe_dofv(n = n, df = df, ncp = 0, seed = seed + 1L)
  sim_ofv <- rep(1000, n)
  alt1_ofv <- sim_ofv + power_d
  alt2_ofv <- sim_ofv - type1_d

  fit <- fake_sse_fit()
  rows <- list()
  for (i in seq_len(n)) {
    rows[[length(rows) + 1L]] <- nlmixr2utils::rawResultsRow(
      fit, source = "sse", hypothesis = "simulation", sample = i,
      modelLabel = "fake_sse_fit", role = "simulation",
      theta = c(tka = 0.55, tcl = 1.05),
      omega = c("omega(eta.ka,eta.ka)" = 0.31),
      objf = sim_ofv[[i]]
    )
    rows[[length(rows) + 1L]] <- nlmixr2utils::rawResultsRow(
      fit, source = "sse", hypothesis = "alternative_1", sample = i,
      modelLabel = "alt1", role = "alternative",
      theta = c(tka = 0.52, tcl = 1.04),
      objf = alt1_ofv[[i]]
    )
    rows[[length(rows) + 1L]] <- nlmixr2utils::rawResultsRow(
      fit, source = "sse", hypothesis = "alternative_2", sample = i,
      modelLabel = "alt2", role = "alternative",
      theta = c(tka = 0.53, tcl = 1.06),
      objf = alt2_ofv[[i]]
    )
  }
  raw_results <- do.call(rbind, rows)

  alt_specs <- list(
    list(label = "alt1", role = "alternative", hypothesis = "alternative_1",
         est = "focei", schema = list(thetaCols = c("tka", "tcl"),
                                       omegaCols = character(0), sigmaCols = character(0))),
    list(label = "alt2", role = "alternative", hypothesis = "alternative_2",
         est = "focei", schema = list(thetaCols = c("tka", "tcl"),
                                       omegaCols = character(0), sigmaCols = character(0)))
  )
  sim_spec <- fake_sse_fit_specs("fake_sse_fit", "alt1")[[1L]]

  run_info <- list(
    fitName = "fake_sse_fit", samples = n, seed = as.integer(seed),
    parameterSource = "fixed", estimateSimulation = TRUE,
    studySampleSize = subjects, studySampleUnit = "subjects",
    studyIdColumn = "ID", studyObservationCount = 2L * subjects,
    fitSpecs = c(list(sim_spec), alt_specs),
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
    alternativeSpecs = list(
      list(label = "alt1", est = "focei", control = list(print = 0L),
           isFit = TRUE, hasDataOverride = FALSE),
      list(label = "alt2", est = "focei", control = list(print = 0L),
           isFit = TRUE, hasDataOverride = FALSE)
    ),
    outputDir = tempdir(), timestamp = Sys.time(),
    referenceValues = outputs$referenceValues,
    initialValues = outputs$initialValues,
    parameterSummary = outputs$parameterSummary,
    ofvSummary = outputs$ofvSummary,
    powerSummary = outputs$powerSummary
  )
}

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
      # binding_nu is a per-OMEGA-BLOCK quantity (the inverse-Wishart degrees
      # of freedom implied by matching this block's reported SEs), not a
      # per-parameter one -- but this fixture's OMEGA block has exactly one
      # eta, so the block-level nu and this row's nu coincide. Computed with
      # the SAME formula .omegaWishartSpec() uses (p + 3 + 2*(mean/sd)^2),
      # here with p = 1. NA for THETA rows, which are Normal-drawn and have no
      # inverse-Wishart degrees of freedom at all.
      binding_nu = c(NA_real_, 1 + 3 + 2 * (omega / se)^2),
      stringsAsFactors = FALSE
    )
  )
  sse
}
