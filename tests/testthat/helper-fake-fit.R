`$.nlmixr2SSEFakeFit` <- function(x, name) {
  x[[name]]
}

# The fake fit must claim `nlmixr2FitCore` to satisfy .validateSSEFit(), which
# means nlmixr2est's `$.nlmixr2FitCore` is a candidate method for it -- and that
# one assumes a real fit environment, so it errors on a plain list. Defining
# `$.nlmixr2SSEFakeFit` here is not enough on its own: package code calling
# `fit$theta` dispatches from inside the nlmixr2sse namespace, where a method
# defined only in the test environment is invisible. nlmixr2est's method then
# wins, errors, and callers that wrap field access in tryCatch (such as
# .rawResultsSchemaForFit()) silently see NULL and build an empty schema.
# Registering it puts the method in the S3 table so dispatch finds it wherever
# the call originates.
registerS3method("$", "nlmixr2SSEFakeFit", `$.nlmixr2SSEFakeFit`)

fake_sse_fit <- function() {
  dummy_ui <- structure(
    list(iniDf = data.frame(stringsAsFactors = FALSE)),
    class = "rxUi"
  )
  fit <- list(
    finalUiEnv = dummy_ui,
    ui = dummy_ui,
    est = "focei",
    control = list(print = 0L),
    theta = c(tka = 0.45, tcl = 1),
    omega = matrix(
      0.3,
      1L,
      1L,
      dimnames = list("eta.ka", "eta.ka")
    ),
    sigma = matrix(numeric(0), 0L, 0L),
    parFixedDf = data.frame(
      Estimate = c(0.45, 1),
      `Back-transformed` = c(exp(0.45), exp(1)),
      SE = c(0.1, 0.2),
      row.names = c("tka", "tcl"),
      check.names = FALSE
    ),
    iniDf = data.frame(
      name = c("tka", "tcl", "eta.ka"),
      ntheta = c(1L, 2L, NA),
      neta1 = c(NA, NA, 1L),
      neta2 = c(NA, NA, 1L),
      lower = c(-Inf, -Inf, NA),
      upper = c(Inf, Inf, NA),
      fix = c(FALSE, FALSE, FALSE),
      stringsAsFactors = FALSE
    ),
    cov = matrix(
      c(0.04, 0.01, 0.01, 0.09),
      2L,
      2L,
      dimnames = list(c("tka", "tcl"), c("tka", "tcl"))
    ),
    objf = 123.45
  )
  class(fit) <- c("nlmixr2SSEFakeFit", "nlmixr2FitCore")
  fit
}

fake_sse_fit_specs <- function(fit_label = "fake_sse_fit", alt_label = "alt1") {
  list(
    list(
      label = fit_label,
      role = "simulation",
      hypothesis = "simulation",
      est = "focei",
      schema = list(
        thetaCols = c("tka", "tcl"),
        omegaCols = "omega(eta.ka,eta.ka)",
        sigmaCols = character(0)
      )
    ),
    list(
      label = alt_label,
      role = "alternative",
      hypothesis = "alternative_1",
      est = "focei",
      schema = list(
        thetaCols = c("tka", "tcl"),
        omegaCols = character(0),
        sigmaCols = character(0)
      )
    )
  )
}

fake_sse_initial_estimates <- function(samples = 2L) {
  data.frame(
    sample = seq_len(samples),
    tka = c(0.50, 0.70)[seq_len(samples)],
    tcl = c(1.00, 1.20)[seq_len(samples)],
    "omega(eta.ka,eta.ka)" = c(0.30, 0.32)[seq_len(samples)],
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

fake_sse_raw_results <- function(
  fit = unclass(fake_sse_fit()),
  fit_label = "fake_sse_fit",
  alt_label = "alt1"
) {
  do.call(
    rbind,
    list(
      nlmixr2utils::rawResultsRow(
        fit,
        source = "sse",
        hypothesis = "simulation",
        sample = 1L,
        modelLabel = fit_label,
        role = "simulation",
        theta = c(tka = 0.55, tcl = 1.05),
        omega = c("omega(eta.ka,eta.ka)" = 0.31),
        objf = 100
      ),
      nlmixr2utils::rawResultsRow(
        fit,
        source = "sse",
        hypothesis = "alternative_1",
        sample = 1L,
        modelLabel = alt_label,
        role = "alternative",
        theta = c(tka = 0.52, tcl = 1.04),
        objf = 103
      ),
      nlmixr2utils::rawResultsRow(
        fit,
        source = "sse",
        hypothesis = "simulation",
        sample = 2L,
        modelLabel = fit_label,
        role = "simulation",
        theta = c(tka = 0.72, tcl = 1.18),
        omega = c("omega(eta.ka,eta.ka)" = 0.33),
        objf = 101
      ),
      nlmixr2utils::rawResultsRow(
        fit,
        source = "sse",
        hypothesis = "alternative_1",
        sample = 2L,
        modelLabel = alt_label,
        role = "alternative",
        theta = c(tka = 0.68, tcl = 1.16),
        objf = 105
      )
    )
  )
}

write_fake_sse_run_dir <- function(
  dir,
  fit_label = "fake_sse_fit",
  alt_label = "alt1",
  completed = TRUE
) {
  fit <- unclass(fake_sse_fit())
  raw_results <- fake_sse_raw_results(
    fit = fit,
    fit_label = fit_label,
    alt_label = alt_label
  )
  initial_wide <- fake_sse_initial_estimates()

  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dir, "simulations"), showWarnings = FALSE)
  dir.create(file.path(dir, "fits"), showWarnings = FALSE)
  nlmixr2utils::writeRawResults(raw_results, dir)
  utils::write.csv(
    initial_wide,
    file.path(dir, "initial_estimates.csv"),
    row.names = FALSE,
    na = "NA"
  )
  saveRDS(
    list(
      runInfo = list(status = if (completed) "completed" else "running"),
      referenceValues = .wideToReferenceValues(initial_wide),
      initialValues = .wideToLongValues(initial_wide),
      parameterSummary = .emptyParameterSummary(),
      ofvSummary = .emptyOfvSummary(),
      powerSummary = .emptyOfvSummary()
    ),
    file.path(dir, "sse_summary.rds")
  )
  saveRDS(
    list(
      status = if (completed) "completed" else "running",
      fitName = fit_label,
      samples = 2L,
      seed = 123L,
      parameterSource = "fixed",
      estimateSimulation = TRUE,
      studySampleSize = 12L,
      studySampleUnit = "subjects",
      studyIdColumn = "ID",
      studyObservationCount = 24L,
      projectedTotalFits = 4L,
      alternativeLabels = alt_label,
      modelLabels = c(fit_label, alt_label),
      fitSpecs = fake_sse_fit_specs(fit_label, alt_label),
      alternativeSpecs = list(list(
        label = alt_label,
        est = "focei",
        control = list(print = 0L),
        isFit = TRUE,
        hasDataOverride = FALSE
      )),
      control = runSSEControl(workers = 1L),
      schemaVersion = 1L,
      startedAt = Sys.time(),
      completedAt = if (completed) Sys.time() else NULL
    ),
    file.path(dir, "run_info.rds")
  )
}

fake_sse_object <- function(
  fit_label = "fake_sse_fit",
  alt_label = "alt1"
) {
  initial_wide <- fake_sse_initial_estimates()
  raw_results <- fake_sse_raw_results(
    fit = unclass(fake_sse_fit()),
    fit_label = fit_label,
    alt_label = alt_label
  )
  outputs <- .computeSSEOutputs(
    rawResults = raw_results,
    initialWide = initial_wide,
    fitSpecsSnapshot = fake_sse_fit_specs(fit_label, alt_label),
    runInfo = list(
      fitName = fit_label,
      samples = 2L,
      parameterSource = "fixed",
      estimateSimulation = TRUE,
      studySampleSize = 12L,
      studySampleUnit = "subjects",
      studyIdColumn = "ID",
      studyObservationCount = 24L,
      fitSpecs = fake_sse_fit_specs(fit_label, alt_label),
      control = runSSEControl(workers = 1L)
    )
  )
  .newNlmixr2SSE(
    runInfo = list(
      fitName = fit_label,
      samples = 2L,
      parameterSource = "fixed",
      estimateSimulation = TRUE,
      studySampleSize = 12L,
      studySampleUnit = "subjects",
      studyIdColumn = "ID",
      studyObservationCount = 24L,
      fitSpecs = fake_sse_fit_specs(fit_label, alt_label),
      control = runSSEControl(workers = 1L)
    ),
    rawResults = raw_results,
    alternativeSpecs = list(list(
      label = alt_label,
      est = "focei",
      control = list(print = 0L),
      isFit = TRUE,
      hasDataOverride = FALSE
    )),
    outputDir = tempdir(),
    timestamp = Sys.time(),
    referenceValues = outputs$referenceValues,
    initialValues = outputs$initialValues,
    parameterSummary = outputs$parameterSummary,
    ofvSummary = outputs$ofvSummary,
    powerSummary = outputs$powerSummary
  )
}

capture_sse_error <- function(expr) {
  tryCatch(
    {
      force(expr)
      NULL
    },
    error = function(e) e
  )
}
