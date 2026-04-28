#' Run fixed-parameter stochastic simulation and estimation
#'
#' `runSSE()` simulates `samples` datasets from a fitted `nlmixr2` reference
#' model, optionally refits the simulation model itself, refits any alternative
#' hypotheses supplied via [sseModel()], writes canonical SSE artifacts, and
#' returns an `nlmixr2SSE` object containing raw results plus tidy parameter and
#' OFV summaries.
#'
#' Phase 6 adds resumable run-state tracking, add-models support on completed
#' SSE directories, and summary recomputation via [recomputeSSE()], while
#' retaining the fixed, raw-results-driven, and covariance-based simulation
#' modes from Phase 5.
#'
#' @param fit Reference fitted `nlmixr2` model.
#' @param alternativeModels Optional `sseModel()` object or list of
#'   `sseModel()` objects.
#' @param samples Integer number of SSE replicates to run.
#' @param seed Optional master seed for the run.
#' @param control SSE control object created by [runSSEControl()].
#' @param outputDir Optional output directory. When `NULL`, a numbered
#'   `<fitName>_sse_<N>` directory is created or resumed.
#' @param restart Logical flag controlling overwrite/resume behavior through the
#'   shared run-directory helper.
#' @param fitName Optional fit label used for the run directory and simulation
#'   model metadata. When `NULL`, the label is derived from the `fit`
#'   expression.
#' @param ... Reserved for future extensions. Unsupported arguments raise an
#'   error.
#'
#' @return An object of class `c("nlmixr2SSE", "list")`.
#' @export
runSSE <- function(
  fit,
  alternativeModels = NULL,
  samples = 100L,
  seed = NULL,
  control = runSSEControl(),
  outputDir = NULL,
  restart = FALSE,
  fitName = NULL,
  ...
) {
  dots <- list(...)
  if (length(dots) > 0L) {
    .abortSSE("Unsupported SSE argument(s): {names(dots)}.")
  }

  .validateSSEFit(fit)
  checkmate::assertIntegerish(
    samples,
    len = 1L,
    any.missing = FALSE,
    lower = 1
  )
  if (!is.null(seed)) {
    checkmate::assertNumber(seed, finite = TRUE)
  }
  checkmate::assertFlag(restart)
  if (!inherits(control, "nlmixr2SSEControl")) {
    .abortSSE("{.arg control} must be created by {.fn runSSEControl}.")
  }
  if (isTRUE(control$addModels) && isTRUE(restart)) {
    .abortSSE(
      "{.arg restart = TRUE} cannot be combined with {.arg control$addModels = TRUE}."
    )
  }
  if (!is.null(outputDir)) {
    checkmate::assertString(outputDir)
  }
  if (!is.null(fitName)) {
    .validateModelLabel(fitName, arg = "fitName")
  } else {
    fitName <- nlmixr2utils::deriveFitName(substitute(fit))
  }

  alternatives <- .normalizeAlternativeModels(alternativeModels, fitName)
  if (
    !isTRUE(control$addModels) &&
      !isTRUE(control$estimateSimulation) &&
      length(alternatives) == 0L
  ) {
    .abortSSE(
      "Supply at least one alternative model when {.arg estimateSimulation = FALSE}."
    )
  }
  if (any(vapply(alternatives, function(x) !is.null(x$data), logical(1)))) {
    cli::cli_inform(c(
      "i" = "Alternative-model {.arg data} overrides are recorded but ignored in the current SSE phase."
    ))
  }

  samples <- as.integer(samples)
  run_dir <- nlmixr2utils::resolveRunDir(
    "sse",
    fitName,
    restart,
    outputDir = outputDir
  )
  if (isTRUE(control$addModels) && identical(run_dir$mode, "new")) {
    .abortSSE(
      "Adding models requires an existing completed SSE run directory."
    )
  }
  if (identical(run_dir$mode, "overwrite") && dir.exists(run_dir$path)) {
    unlink(run_dir$path, recursive = TRUE, force = TRUE)
  }
  dir.create(run_dir$path, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(run_dir$path, "simulations"), showWarnings = FALSE)
  dir.create(file.path(run_dir$path, "fits"), showWarnings = FALSE)

  abs_output_dir <- normalizePath(run_dir$path, mustWork = FALSE)
  master_seed <- nlmixr2utils::withRunSeed(
    abs_output_dir,
    seed = seed,
    prefix = "sse"
  )

  existing_run_info <- .readRunInfo(abs_output_dir)
  existing_state <- .readSSEState(abs_output_dir)
  data_cache <- nlmixr2utils::taskCache(abs_output_dir, "sse_data")
  fit_cache <- nlmixr2utils::taskCache(abs_output_dir, "sse_fit_rows")
  existing_state <- .syncSSEState(
    dir = abs_output_dir,
    state = existing_state,
    dataCache = data_cache,
    fitCache = fit_cache,
    status = existing_state$status,
    phase = existing_state$phase,
    fitName = if (is.na(existing_state$fitName)) {
      fitName
    } else {
      existing_state$fitName
    },
    samples = if (is.na(existing_state$samples)) {
      samples
    } else {
      existing_state$samples
    },
    modelLabels = existing_state$modelLabels,
    fitFailures = existing_state$fitFailures,
    addRuns = existing_state$addRuns
  )
  existing_raw_results <- if (
    file.exists(file.path(abs_output_dir, "raw_results.rds")) ||
      file.exists(file.path(abs_output_dir, "raw_results.csv"))
  ) {
    nlmixr2utils::readRawResults(abs_output_dir)
  } else {
    NULL
  }

  sim_schema <- .rawResultsSchemaForFit(fit)
  worker_desc <- .workerDescription(control$workers)

  if (isTRUE(control$addModels)) {
    if (
      !is.list(existing_run_info) ||
        !identical(existing_run_info$status, "completed")
    ) {
      .abortSSE(
        "{.arg control$addModels = TRUE} requires an existing completed SSE run."
      )
    }
    if (!identical(as.integer(existing_run_info$samples), samples)) {
      .abortSSE(
        "Add-models requests must use the original {.arg samples = {existing_run_info$samples}}."
      )
    }
    if (
      !identical(existing_run_info$parameterSource, control$parameterSource)
    ) {
      .abortSSE(
        "Add-models requests must reuse the original {.arg parameterSource = {existing_run_info$parameterSource}}."
      )
    }
    old_control <- existing_run_info$control %||% list()
    if (
      !identical(
        isTRUE(old_control$randomEstimationInits),
        isTRUE(control$randomEstimationInits)
      )
    ) {
      .abortSSE(
        "Add-models requests must reuse the original {.arg randomEstimationInits} setting."
      )
    }
    if (!identical(isTRUE(old_control$updateFix), isTRUE(control$updateFix))) {
      .abortSSE(
        "Add-models requests must reuse the original {.arg updateFix} setting."
      )
    }
    alternatives <- .selectAddModelAlternatives(alternatives, existing_run_info)
  } else if (!is.null(existing_run_info)) {
    .validateResumeRequest(
      existingRunInfo = existing_run_info,
      fitName = fitName,
      samples = samples,
      control = control,
      requestedLabels = .requestedModelLabels(.fitTaskSpecs(
        fit,
        fitName,
        alternatives,
        control
      ))
    )
    if (
      identical(existing_run_info$status, "completed") &&
        .completedSSEArtifactsExist(abs_output_dir)
    ) {
      return(.loadCompletedSSE(
        outputDir = abs_output_dir,
        runInfo = existing_run_info
      ))
    }
  }

  task_control <- control
  if (isTRUE(control$addModels)) {
    task_control$estimateSimulation <- FALSE
  }
  fit_specs <- .fitTaskSpecs(fit, fitName, alternatives, task_control)
  fit_specs_snapshot <- .fitSpecsSnapshot(fit_specs)
  existing_fit_specs <- if (isTRUE(control$addModels)) {
    .fitSpecsFromRunInfo(existing_run_info, rawResults = existing_raw_results)
  } else {
    list()
  }
  all_fit_specs <- if (isTRUE(control$addModels)) {
    c(existing_fit_specs, fit_specs_snapshot)
  } else {
    fit_specs_snapshot
  }
  projected_total <- samples * length(all_fit_specs)
  pending_projected <- samples * length(fit_specs)

  schema_inputs <- c(
    if (!is.null(existing_raw_results)) {
      list(.rawResultsSchemaFromHeader(attr(
        existing_raw_results,
        "rawResultsHeader",
        exact = TRUE
      )))
    } else {
      list(sim_schema)
    },
    lapply(fit_specs, `[[`, "schema")
  )
  union_schema <- .mergeRawResultsSchemas(schema_inputs)

  cli::cli_rule(left = "SSE")
  cli::cli_inform(c(
    "i" = "Model             : {fitName}",
    "i" = "Parameter source  : {control$parameterSource}",
    "i" = "Replicates        : {samples}",
    "i" = "Projected fits    : {projected_total}",
    "i" = "Scheduled fits    : {pending_projected}",
    "i" = "Execution         : {worker_desc}",
    "i" = "Run directory     : {abs_output_dir}"
  ))
  if (projected_total > 500L) {
    cli::cli_warn(c(
      "!" = "This SSE run projects {projected_total} total fits.",
      "i" = "Adjust {.arg samples} or {.arg alternativeModels} if you want a smaller run."
    ))
  }
  if (identical(control$parameterSource, "covariance")) {
    cli::cli_inform(c(
      "i" = "Thetas are drawn from {.field fit$cov}; OMEGA and SIGMA stay at the fitted point estimates."
    ))
  }
  cli::cli_rule()

  ref_data <- .getFitData(fit)
  study_info <- .inferStudySampleInfo(ref_data)

  run_info <- existing_run_info %||% list()
  run_info$status <- "running"
  run_info$fitName <- fitName
  run_info$samples <- samples
  run_info$seed <- master_seed
  run_info$parameterSource <- if (isTRUE(control$addModels)) {
    existing_run_info$parameterSource
  } else {
    control$parameterSource
  }
  run_info$estimateSimulation <- if (isTRUE(control$addModels)) {
    isTRUE(existing_run_info$estimateSimulation)
  } else {
    control$estimateSimulation
  }
  run_info$projectedTotalFits <- projected_total
  run_info$workers <- control$workers
  run_info$runDirMode <- run_dir$mode
  run_info$alternativeLabels <- vapply(
    Filter(function(spec) identical(spec$role, "alternative"), all_fit_specs),
    `[[`,
    character(1),
    "label"
  )
  run_info$modelLabels <- vapply(all_fit_specs, `[[`, character(1), "label")
  run_info$fitSpecs <- all_fit_specs
  run_info$alternativeSpecs <- if (isTRUE(control$addModels)) {
    c(
      .alternativeSnapshotsFromRunInfo(existing_run_info %||% list()),
      .alternativeSnapshots(alternatives)
    )
  } else {
    .alternativeSnapshots(alternatives)
  }
  run_info$schemaVersion <- union_schema$schemaVersion
  run_info$parameterSourceInfo <- if (isTRUE(control$addModels)) {
    existing_run_info$parameterSourceInfo
  } else {
    .resolveSimulationParameters(
      fit = fit,
      samples = samples,
      control = control,
      outputDir = abs_output_dir,
      schema = sim_schema
    )$info
  }
  run_info$studySampleSize <- if (isTRUE(control$addModels)) {
    existing_run_info$studySampleSize %||% study_info$studySampleSize
  } else {
    study_info$studySampleSize
  }
  run_info$studySampleUnit <- if (isTRUE(control$addModels)) {
    existing_run_info$studySampleUnit %||% study_info$studySampleUnit
  } else {
    study_info$studySampleUnit
  }
  run_info$studyIdColumn <- if (isTRUE(control$addModels)) {
    existing_run_info$studyIdColumn %||% study_info$studyIdColumn
  } else {
    study_info$studyIdColumn
  }
  run_info$studyObservationCount <- if (isTRUE(control$addModels)) {
    existing_run_info$studyObservationCount %||% study_info$studyObservationCount
  } else {
    study_info$studyObservationCount
  }
  run_info$control <- control
  run_info$control$estimateSimulation <- run_info$estimateSimulation
  run_info$control$parameterSource <- run_info$parameterSource
  run_info$call <- match.call()
  run_info$startedAt <- run_info$startedAt %||% Sys.time()
  saveRDS(run_info, file.path(abs_output_dir, "run_info.rds"))
  state <- .syncSSEState(
    dir = abs_output_dir,
    state = existing_state,
    dataCache = data_cache,
    fitCache = fit_cache,
    status = "running",
    phase = if (isTRUE(control$addModels)) "add-models" else "simulation",
    fitName = fitName,
    samples = samples,
    modelLabels = run_info$modelLabels,
    fitFailures = NA_integer_,
    addRuns = existing_state$addRuns
  )

  sample_keys <- as.character(seq_len(samples))
  if (isTRUE(control$addModels)) {
    initial_wide <- .readInitialEstimates(abs_output_dir)
    sim_records <- lapply(
      seq_len(samples),
      function(sample_id) {
        .readSavedSimulationRecord(
          outputDir = abs_output_dir,
          sample = sample_id,
          simSchema = sim_schema,
          fit = fit,
          initialWide = initial_wide,
          dataCache = data_cache
        )
      }
    )
    for (record in sim_records) {
      key <- as.character(record$sample)
      if (!data_cache$has(key)) {
        data_cache$put(key, record)
      }
    }
    state <- .syncSSEState(
      dir = abs_output_dir,
      state = state,
      dataCache = data_cache,
      fitCache = fit_cache,
      status = "running",
      phase = "fitting",
      fitName = fitName,
      samples = samples,
      modelLabels = run_info$modelLabels,
      fitFailures = NA_integer_,
      addRuns = state$addRuns
    )
    pending_data <- integer(0)
  } else {
    parameter_source <- .resolveSimulationParameters(
      fit = fit,
      samples = samples,
      control = control,
      outputDir = abs_output_dir,
      schema = sim_schema
    )
    sim_param_sets <- parameter_source$records
    pending_data <- as.integer(
      nlmixr2utils::pendingTasks(sample_keys, data_cache$keys())
    )

    if (length(pending_data) > 0L) {
      nlmixr2utils::.withWorkerPlan(control$workers, {
        nlmixr2utils::.plap(
          pending_data,
          function(sample_id) {
            record <- .simulationRecord(
              sample = sample_id,
              fit = fit,
              refData = ref_data,
              control = control,
              outputDir = abs_output_dir,
              schema = sim_schema,
              paramSet = sim_param_sets[[sample_id]]
            )
            data_cache$put(as.character(sample_id), record)
            invisible(NULL)
          },
          .label = function(sample_id) {
            paste0("simulation ", sample_id, "/", samples)
          }
        )
      })
    }

    state <- .syncSSEState(
      dir = abs_output_dir,
      state = state,
      dataCache = data_cache,
      fitCache = fit_cache,
      status = "running",
      phase = "fitting",
      fitName = fitName,
      samples = samples,
      modelLabels = run_info$modelLabels,
      fitFailures = NA_integer_,
      addRuns = state$addRuns
    )
    sim_records <- lapply(sample_keys, data_cache$get)
    initial_wide <- .initialEstimatesWide(
      sim_records,
      sim_schema,
      seq_len(samples)
    )
    .writeInitialEstimates(initial_wide, abs_output_dir)
  }

  fit_keys <- unlist(
    lapply(fit_specs, function(spec) {
      vapply(
        seq_len(samples),
        function(i) .fitTaskKey(spec$label, i),
        character(1)
      )
    }),
    use.names = FALSE
  )
  pending_fit <- nlmixr2utils::pendingTasks(fit_keys, fit_cache$keys())

  if (length(pending_fit) > 0L) {
    label_index <- stats::setNames(
      seq_along(fit_specs),
      vapply(
        fit_specs,
        `[[`,
        character(1),
        "label"
      )
    )

    nlmixr2utils::.withWorkerPlan(control$workers, {
      nlmixr2utils::.plap(
        pending_fit,
        function(key) {
          parts <- strsplit(key, "::", fixed = TRUE)[[1L]]
          label <- parts[[1L]]
          sample_id <- as.integer(parts[[2L]])
          spec <- fit_specs[[label_index[[label]]]]
          sim_record <- data_cache$get(as.character(sample_id))
          fit_record <- .fitTaskRecord(
            sample = sample_id,
            spec = spec,
            simRecord = sim_record,
            unionSchema = union_schema,
            outputDir = abs_output_dir,
            saveFits = control$saveFits,
            control = control
          )
          fit_cache$put(key, fit_record)
          invisible(NULL)
        },
        .label = function(key) {
          parts <- strsplit(key, "::", fixed = TRUE)[[1L]]
          paste0(parts[[1L]], " ", parts[[2L]], "/", samples)
        }
      )
    })
  }

  fit_records <- lapply(fit_keys, fit_cache$get)
  new_rows <- do.call(rbind, lapply(fit_records, `[[`, "row"))
  raw_results <- if (
    isTRUE(control$addModels) && !is.null(existing_raw_results)
  ) {
    rbind(existing_raw_results, new_rows)
  } else {
    new_rows
  }
  raw_results <- raw_results[
    order(
      raw_results$sample,
      raw_results$role,
      raw_results$model_label
    ),
    ,
    drop = FALSE
  ]

  raw_results_out <- nlmixr2utils::writeRawResults(raw_results, abs_output_dir)
  raw_results <- raw_results_out$data
  if (isTRUE(control$addModels)) {
    add_idx <- .nextArtifactIndex(abs_output_dir, "raw_results_add", "csv")
    add_basename <- paste0("raw_results_add", add_idx)
    add_raw_results_out <- nlmixr2utils::writeRawResults(
      new_rows,
      abs_output_dir,
      basename = add_basename
    )
    .writeSseResults(
      .computeParameterSummary(
        rawResults = add_raw_results_out$data,
        initialWide = initial_wide,
        specMap = .specMapFromSnapshots(fit_specs_snapshot),
        outFilter = control$outFilter
      ),
      abs_output_dir,
      basename = paste0("sse_results_add", add_idx)
    )
    state$addRuns <- c(
      state$addRuns,
      list(list(
        index = add_idx,
        labels = vapply(fit_specs_snapshot, `[[`, character(1), "label"),
        rawResultsBasename = add_basename,
        resultsBasename = paste0("sse_results_add", add_idx),
        completedAt = Sys.time()
      ))
    )
  }

  outputs <- .computeSSEOutputs(
    rawResults = raw_results,
    initialWide = initial_wide,
    fitSpecsSnapshot = all_fit_specs,
    runInfo = run_info,
    outFilter = control$outFilter,
    refOfv = control$refOfv
  )

  n_failed <- sum(!is.na(raw_results$error_message))
  completed_info <- run_info
  completed_info$status <- "completed"
  completed_info$completedAt <- Sys.time()
  completed_info$rawResultsRows <- nrow(raw_results)
  completed_info$fitFailures <- as.integer(n_failed)
  completed_info$addRuns <- state$addRuns

  .writeSseResults(outputs$parameterSummary, abs_output_dir)
  .writeSseSummary(
    runInfo = completed_info,
    referenceValues = outputs$referenceValues,
    initialValues = outputs$initialValues,
    parameterSummary = outputs$parameterSummary,
    ofvSummary = outputs$ofvSummary,
    dir = abs_output_dir
  )

  saveRDS(completed_info, file.path(abs_output_dir, "run_info.rds"))
  .syncSSEState(
    dir = abs_output_dir,
    state = state,
    dataCache = data_cache,
    fitCache = fit_cache,
    status = "completed",
    phase = "completed",
    fitName = fitName,
    samples = samples,
    modelLabels = run_info$modelLabels,
    fitFailures = as.integer(n_failed),
    addRuns = state$addRuns
  )

  cli::cli_rule(left = "SSE complete")
  if (n_failed == 0L) {
    cli::cli_inform(c(
      "v" = "All {projected_total} fit(s) completed without caught exceptions.",
      "i" = "Results saved to: {abs_output_dir}"
    ))
  } else {
    cli::cli_warn(c(
      "!" = "{n_failed} fit(s) raised caught exceptions during the SSE run.",
      "i" = "Their raw-results rows carry {.field error_message} text and saved fit artifacts use class {.cls nlmixr2SSEFitError}."
    ))
    cli::cli_inform(c(
      "i" = "Results saved to: {abs_output_dir}"
    ))
  }
  cli::cli_rule()

  .newNlmixr2SSE(
    runInfo = completed_info,
    rawResults = raw_results,
    alternativeSpecs = completed_info$alternativeSpecs,
    outputDir = abs_output_dir,
    timestamp = completed_info$completedAt,
    referenceValues = outputs$referenceValues,
    initialValues = outputs$initialValues,
    parameterSummary = outputs$parameterSummary,
    ofvSummary = outputs$ofvSummary,
    powerSummary = outputs$powerSummary
  )
}
