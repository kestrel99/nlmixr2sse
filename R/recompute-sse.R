# Phase 6 restart and recompute helpers ---------------------------------------

.sseStatePrototype <- function() {
  list(
    status = "new",
    phase = "new",
    fitName = NA_character_,
    samples = NA_integer_,
    modelLabels = character(0),
    completedSamples = integer(0),
    completedFitTasks = character(0),
    fitFailures = NA_integer_,
    addRuns = list(),
    updatedAt = Sys.time()
  )
}

.readSSEState <- function(dir) {
  state <- nlmixr2utils::readRunState(dir, "sse")
  out <- .sseStatePrototype()
  if (is.null(state)) {
    return(out)
  }

  for (nm in intersect(names(out), names(state))) {
    out[[nm]] <- state[[nm]]
  }
  out$completedSamples <- sort(unique(as.integer(out$completedSamples)))
  out$completedFitTasks <- unique(as.character(out$completedFitTasks))
  out$modelLabels <- unname(as.character(out$modelLabels))
  out$addRuns <- out$addRuns %||% list()
  out
}

.writeSSEState <- function(dir, state) {
  state$completedSamples <- sort(unique(as.integer(state$completedSamples)))
  state$completedFitTasks <- unique(as.character(state$completedFitTasks))
  state$modelLabels <- unname(as.character(state$modelLabels))
  state$updatedAt <- Sys.time()
  nlmixr2utils::writeRunState(dir, state = state, schema = "sse")
  invisible(state)
}

.cacheSamples <- function(keys) {
  keys <- suppressWarnings(as.integer(keys))
  sort(unique(keys[!is.na(keys)]))
}

.syncSSEState <- function(
  dir,
  state = .readSSEState(dir),
  dataCache = NULL,
  fitCache = NULL,
  status = state$status,
  phase = state$phase,
  fitName = state$fitName,
  samples = state$samples,
  modelLabels = state$modelLabels,
  fitFailures = state$fitFailures,
  addRuns = state$addRuns
) {
  if (!is.null(dataCache)) {
    state$completedSamples <- .cacheSamples(dataCache$keys())
  }
  if (!is.null(fitCache)) {
    state$completedFitTasks <- sort(unique(as.character(fitCache$keys())))
  }
  state$status <- status
  state$phase <- phase
  state$fitName <- fitName
  state$samples <- as.integer(samples)
  state$modelLabels <- unname(as.character(modelLabels))
  state$fitFailures <- as.integer(fitFailures)
  state$addRuns <- addRuns %||% list()
  .writeSSEState(dir, state)
}

.readRunInfo <- function(dir) {
  path <- file.path(dir, "run_info.rds")
  if (!file.exists(path)) {
    return(NULL)
  }
  readRDS(path)
}

.completedSSEArtifactsExist <- function(dir) {
  all(file.exists(file.path(
    dir,
    c("raw_results.rds", "initial_estimates.csv", "sse_summary.rds")
  )))
}

.requestedModelLabels <- function(fitSpecs) {
  vapply(fitSpecs, `[[`, character(1), "label")
}

.existingModelLabels <- function(runInfo, rawResults = NULL) {
  labels <- runInfo$modelLabels %||% character(0)
  labels <- unname(as.character(labels))
  if (length(labels) > 0L) {
    return(labels)
  }
  if (is.null(rawResults) || nrow(rawResults) == 0L) {
    return(character(0))
  }
  unique(rawResults$model_label)
}

.validateResumeRequest <- function(
  existingRunInfo,
  fitName,
  samples,
  control,
  requestedLabels,
  rxThreads = NULL
) {
  if (is.null(existingRunInfo)) {
    return(invisible(NULL))
  }

  if (!identical(existingRunInfo$fitName, fitName)) {
    .abortSSE(
      "Existing run directory belongs to {.val {existingRunInfo$fitName}}, not {.val {fitName}}."
    )
  }
  if (!identical(as.integer(existingRunInfo$samples), as.integer(samples))) {
    .abortSSE(
      "Existing run directory was created with {.arg samples = {existingRunInfo$samples}}, not {.arg samples = {samples}}."
    )
  }
  if (!identical(existingRunInfo$parameterSource, control$parameterSource)) {
    .abortSSE(
      "Existing run directory uses {.arg parameterSource = {existingRunInfo$parameterSource}}, not {.arg parameterSource = {control$parameterSource}}."
    )
  }
  if (
    !identical(
      isTRUE(existingRunInfo$estimateSimulation),
      isTRUE(control$estimateSimulation)
    )
  ) {
    .abortSSE(
      "Existing run directory was created with a different {.arg estimateSimulation} setting."
    )
  }

  if (!is.null(rxThreads)) {
    recorded <- existingRunInfo$rxThreads
    if (is.null(recorded)) {
      cli::cli_warn(c(
        "!" = "This run directory predates {.arg rxThreads} tracking, so the thread count it used is unknown.",
        "i" = "rxode2 thread count changes simulated values, so resumed replicates may not be comparable to the existing ones."
      ))
    } else if (!identical(as.integer(recorded), as.integer(rxThreads))) {
      .abortSSE(
        paste0(
          "Existing run directory used {.arg rxThreads = {as.integer(recorded)}}, ",
          "but this run resolves to {.arg rxThreads = {as.integer(rxThreads)}}. ",
          "rxode2 thread count changes simulated values, so the replicates would not be comparable. ",
          "Use {.code runSSEControl(rxThreads = {as.integer(recorded)})} to resume, or {.code restart = TRUE} for a fresh run."
        )
      )
    }
  }

  if (identical(control$parameterSource, "covariance")) {
    recordedDraw <- existingRunInfo$covarianceDraw
    # NULL and NA both mean "legacy, mode unknown". NA specifically arises after
    # an addModels run against a legacy directory, which records NA by design --
    # without treating it as legacy here, that run would fall into the generic
    # mismatch branch below and report the original mode as "NA".
    if (is.null(recordedDraw) || is.na(recordedDraw)) {
      # A covariance run with no recorded mode necessarily predates OMEGA
      # drawing, so its existing replicates hold OMEGA at the fitted values.
      # Appending OMEGA-varying replicates would put two simulation
      # distributions in one study. This is NOT the rxThreads case, where the
      # thread count is merely unknown but the draw is the same kind -- here the
      # old behaviour is known and known to be incompatible. Abort.
      .abortSSE(
        paste0(
          "This run directory predates {.arg covarianceDraw} tracking, so its ",
          "existing replicates were simulated with OMEGA held at the fitted ",
          "values. Resuming would append replicates that vary OMEGA, mixing two ",
          "simulation distributions in one study. Use {.code restart = TRUE} to ",
          "start a fresh run."
        )
      )
    } else if (!identical(recordedDraw, control$covarianceDraw)) {
      .abortSSE(
        paste0(
          "Existing run directory used {.arg covarianceDraw = {recordedDraw}}, ",
          "but this run requests {.arg covarianceDraw = {control$covarianceDraw}}. ",
          "The draw mode changes every simulated parameter set, so the replicates would not be comparable. ",
          "Use {.code runSSEControl(covarianceDraw = \"{recordedDraw}\")} to resume, or {.code restart = TRUE} for a fresh run."
        )
      )
    }
  }

  existing_labels <- .existingModelLabels(existingRunInfo)
  if (
    length(existing_labels) > 0L &&
      !setequal(existing_labels, unname(as.character(requestedLabels)))
  ) {
    .abortSSE(
      paste0(
        "Requested model set does not match the existing run directory. ",
        "Use `runSSEControl(addModels = TRUE)` to extend a completed run or `restart = TRUE` for a fresh run."
      )
    )
  }

  invisible(NULL)
}

.nextArtifactIndex <- function(dir, stem, extension) {
  pattern <- paste0("^", stem, "([0-9]+)\\.", extension, "$")
  files <- list.files(dir, pattern = pattern, full.names = FALSE)
  if (length(files) == 0L) {
    return(1L)
  }
  idx <- as.integer(sub(pattern, "\\1", files))
  max(idx, na.rm = TRUE) + 1L
}

.savedDatasetPath <- function(dir, sample) {
  file.path(dir, "simulations", sprintf("sim_%04d.rds", sample))
}

.matrixInfoFromTemplate <- function(mat, prefix) {
  if (
    !is.matrix(mat) || length(mat) == 0L || nrow(mat) == 0L || ncol(mat) == 0L
  ) {
    return(data.frame(
      colName = character(0),
      row = integer(0),
      col = integer(0),
      stringsAsFactors = FALSE
    ))
  }

  idx <- which(lower.tri(mat, diag = TRUE), arr.ind = TRUE)
  idx <- idx[order(idx[, "col"], idx[, "row"]), , drop = FALSE]
  row_names <- rownames(mat)
  col_names <- colnames(mat)

  label_one <- function(r, c_idx) {
    row_lab <- if (!is.null(row_names) && length(row_names) >= r) {
      row_names[[r]]
    } else {
      as.character(r)
    }
    col_lab <- if (!is.null(col_names) && length(col_names) >= c_idx) {
      col_names[[c_idx]]
    } else {
      as.character(c_idx)
    }
    paste0(prefix, "(", row_lab, ",", col_lab, ")")
  }

  data.frame(
    colName = vapply(
      seq_len(nrow(idx)),
      function(i) label_one(idx[i, "row"], idx[i, "col"]),
      character(1)
    ),
    row = idx[, "row"],
    col = idx[, "col"],
    stringsAsFactors = FALSE
  )
}

.paramSetFromInitialRow <- function(
  initialWide,
  sample,
  simSchema,
  fit
) {
  row <- initialWide[match(sample, initialWide$sample), , drop = FALSE]
  if (nrow(row) != 1L) {
    .abortSSE(
      "Could not recover initial estimates for saved simulation sample {.val {sample}}."
    )
  }

  theta <- .thetaVectorFromFit(fit)
  if (length(simSchema$thetaCols) > 0L) {
    theta[simSchema$thetaCols] <- as.numeric(row[
      1,
      simSchema$thetaCols,
      drop = TRUE
    ])
  }

  omega <- .fitField(fit, "omega")
  if (!is.matrix(omega)) {
    omega <- matrix(numeric(0), 0L, 0L)
  }
  omega_info <- .matrixInfoFromTemplate(omega, "omega")
  if (nrow(omega_info) > 0L) {
    for (i in seq_len(nrow(omega_info))) {
      value <- as.numeric(row[[omega_info$colName[[i]]]])
      omega[omega_info$row[[i]], omega_info$col[[i]]] <- value
      omega[omega_info$col[[i]], omega_info$row[[i]]] <- value
    }
  }

  sigma <- .fitField(fit, "sigma")
  if (!is.matrix(sigma)) {
    sigma <- matrix(numeric(0), 0L, 0L)
  }
  sigma_info <- .matrixInfoFromTemplate(sigma, "sigma")
  if (nrow(sigma_info) > 0L) {
    for (i in seq_len(nrow(sigma_info))) {
      value <- as.numeric(row[[sigma_info$colName[[i]]]])
      sigma[sigma_info$row[[i]], sigma_info$col[[i]]] <- value
      sigma[sigma_info$col[[i]], sigma_info$row[[i]]] <- value
    }
  }

  list(theta = theta, omega = omega, sigma = sigma)
}

.readSavedSimulationRecord <- function(
  outputDir,
  sample,
  simSchema,
  fit,
  initialWide,
  dataCache = NULL
) {
  path <- .savedDatasetPath(outputDir, sample)
  if (!file.exists(path)) {
    .abortSSE(
      paste0(
        "Saved simulation dataset `",
        basename(path),
        "` was not found. ",
        "`runSSEControl(addModels = TRUE)` requires the original run to keep simulation datasets."
      )
    )
  }

  cached <- if (!is.null(dataCache) && dataCache$has(as.character(sample))) {
    dataCache$get(as.character(sample))
  } else {
    NULL
  }
  param_set <- if (!is.null(cached$paramSet)) {
    cached$paramSet
  } else {
    .paramSetFromInitialRow(initialWide, sample, simSchema, fit)
  }
  initial <- if (!is.null(cached$initial)) {
    cached$initial
  } else {
    as.numeric(initialWide[
      match(sample, initialWide$sample),
      .schemaParamCols(simSchema),
      drop = TRUE
    ])
  }
  names(initial) <- .schemaParamCols(simSchema)

  list(
    sample = as.integer(sample),
    data = readRDS(path),
    initial = initial,
    paramSet = param_set
  )
}

.selectAddModelAlternatives <- function(alternatives, runInfo) {
  if (length(alternatives) == 0L) {
    .abortSSE(
      "Supply at least one new alternative model when {.arg control$addModels = TRUE}."
    )
  }

  existing_labels <- .existingModelLabels(runInfo)
  colliding <- intersect(
    existing_labels,
    vapply(alternatives, `[[`, character(1), "label")
  )
  if (length(colliding) > 0L) {
    .abortSSE(
      "Alternative-model label{?s} {.val {colliding}} already exist in the completed run directory."
    )
  }

  alternatives
}

.computeSSEOutputs <- function(
  rawResults,
  initialWide,
  fitSpecsSnapshot,
  runInfo,
  outFilter = NULL,
  refOfv = NULL
) {
  spec_map <- .specMapFromSnapshots(fitSpecsSnapshot)
  parameter_summary <- .computeParameterSummary(
    rawResults = rawResults,
    initialWide = initialWide,
    specMap = spec_map,
    outFilter = outFilter
  )

  control <- runInfo$control %||% list()
  control$outFilter <- outFilter %||% control$outFilter
  if (!is.null(refOfv)) {
    control$refOfv <- refOfv
  }

  ofv_summary <- .computeOfvSummary(
    rawResults = rawResults,
    fitSpecs = fitSpecsSnapshot,
    control = control,
    initialWide = initialWide
  )

  list(
    referenceValues = .wideToReferenceValues(initialWide),
    initialValues = .wideToLongValues(initialWide),
    parameterSummary = parameter_summary,
    ofvSummary = ofv_summary,
    powerSummary = ofv_summary[!is.na(ofv_summary$direction) & ofv_summary$direction == "power", , drop = FALSE]
  )
}

.writeRecomputeArtifacts <- function(
  dir,
  outputs,
  rawResults,
  runInfo,
  overwrite
) {
  if (isTRUE(overwrite)) {
    .writeSseResults(outputs$parameterSummary, dir, basename = "sse_results")
    .writeSseSummary(
      runInfo = runInfo,
      referenceValues = outputs$referenceValues,
      initialValues = outputs$initialValues,
      parameterSummary = outputs$parameterSummary,
      ofvSummary = outputs$ofvSummary,
      dir = dir,
      basename = "sse_summary"
    )
    return(list(
      resultsBasename = "sse_results",
      summaryBasename = "sse_summary",
      rawResults = rawResults
    ))
  }

  idx <- .nextArtifactIndex(dir, "sse_results_recompute", "csv")
  results_basename <- paste0("sse_results_recompute", idx)
  summary_basename <- paste0("sse_summary_recompute", idx)
  .writeSseResults(outputs$parameterSummary, dir, basename = results_basename)
  .writeSseSummary(
    runInfo = runInfo,
    referenceValues = outputs$referenceValues,
    initialValues = outputs$initialValues,
    parameterSummary = outputs$parameterSummary,
    ofvSummary = outputs$ofvSummary,
    dir = dir,
    basename = summary_basename
  )
  list(
    resultsBasename = results_basename,
    summaryBasename = summary_basename,
    rawResults = rawResults
  )
}

#' Recompute SSE summaries from an existing run directory
#'
#' `recomputeSSE()` rebuilds the parameter and OFV summary artifacts from an
#' existing completed SSE run without rerunning any simulations or fits. This
#' is useful when you want to apply a different `outFilter`, use a different
#' `refOfv`, or recompute the summaries on a subset of the completed
#' replicates.
#'
#' @param runDir Completed SSE run directory.
#' @param outFilter Optional raw-results filter applied before summary
#'   calculations. This accepts the same forms as
#'   [nlmixr2utils::setupRawResultsFilter()].
#' @param refOfv Optional OFV override used when recomputing the OFV summary.
#' @param samples Optional integer number of replicates to include from the
#'   start of the saved run.
#' @param overwrite Logical. When `FALSE`, writes
#'   `sse_results_recompute<k>.csv` and `sse_summary_recompute<k>.rds`.
#'   When `TRUE`, replaces `sse_results.csv` and `sse_summary.rds`.
#'
#' @return An object of class `c("nlmixr2SSE", "list")`.
#' @export
recomputeSSE <- function(
  runDir,
  outFilter = NULL,
  refOfv = NULL,
  samples = NULL,
  overwrite = FALSE
) {
  checkmate::assertString(runDir)
  if (!dir.exists(runDir)) {
    .abortSSE("{.arg runDir} does not exist: {.path {runDir}}.")
  }
  if (!is.null(outFilter)) {
    .validateFilterInput(outFilter, "outFilter")
  }
  if (!is.null(refOfv)) {
    checkmate::assertNumber(refOfv, finite = TRUE)
  }
  if (!is.null(samples)) {
    checkmate::assertIntegerish(
      samples,
      len = 1L,
      any.missing = FALSE,
      lower = 1
    )
    samples <- as.integer(samples)
  }
  checkmate::assertFlag(overwrite)

  abs_output_dir <- normalizePath(runDir, mustWork = FALSE)
  run_info <- .readRunInfo(abs_output_dir)
  if (!is.list(run_info) || !identical(run_info$status, "completed")) {
    .abortSSE(
      "Recomputation requires a completed SSE run with {.file run_info.rds}."
    )
  }

  raw_results <- nlmixr2utils::readRawResults(abs_output_dir)
  initial_wide <- .readInitialEstimates(abs_output_dir)
  if (!is.null(samples)) {
    if (samples > nrow(initial_wide)) {
      .abortSSE(
        "Requested {.arg samples = {samples}} exceeds the {nrow(initial_wide)} saved replicate(s)."
      )
    }
    initial_wide <- initial_wide[initial_wide$sample <= samples, , drop = FALSE]
    raw_results <- raw_results[raw_results$sample <= samples, , drop = FALSE]
    run_info$samples <- samples
  }

  fit_specs <- .fitSpecsFromRunInfo(run_info, rawResults = raw_results)
  outputs <- .computeSSEOutputs(
    rawResults = raw_results,
    initialWide = initial_wide,
    fitSpecsSnapshot = fit_specs,
    runInfo = run_info,
    outFilter = outFilter,
    refOfv = refOfv
  )

  .writeRecomputeArtifacts(
    dir = abs_output_dir,
    outputs = outputs,
    rawResults = raw_results,
    runInfo = run_info,
    overwrite = overwrite
  )

  .newNlmixr2SSE(
    runInfo = run_info,
    rawResults = raw_results,
    alternativeSpecs = .alternativeSnapshotsFromRunInfo(run_info),
    outputDir = abs_output_dir,
    timestamp = Sys.time(),
    referenceValues = outputs$referenceValues,
    initialValues = outputs$initialValues,
    parameterSummary = outputs$parameterSummary,
    ofvSummary = outputs$ofvSummary,
    powerSummary = outputs$powerSummary
  )
}
