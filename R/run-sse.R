#' Initialize an `nlmixr2` stochastic simulation and estimation run
#'
#' `runSSE()` is the main entry point for `nlmixr2sse`. Phase 3 establishes the
#' shared control flow: it validates the reference fit and alternative-model
#' specifications, resolves the numbered run directory via
#' [nlmixr2utils::resolveRunDir()], persists the run seed via
#' [nlmixr2utils::withRunSeed()], and returns an initialized `nlmixr2SSE`
#' object. The fixed-parameter simulation and refit core is added in Phase 4.
#'
#' @param fit Reference fitted `nlmixr2` model.
#' @param alternativeModels Optional `sseModel()` object or list of
#'   `sseModel()` objects.
#' @param samples Integer number of SSE replicates to prepare.
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
  if (!is.null(outputDir)) {
    checkmate::assertString(outputDir)
  }
  if (!is.null(fitName)) {
    .validateModelLabel(fitName, arg = "fitName")
  } else {
    fitName <- nlmixr2utils::deriveFitName(substitute(fit))
  }

  if (identical(control$parameterSource, "covariance")) {
    .validateCovarianceFit(fit)
  }

  alternatives <- .normalizeAlternativeModels(alternativeModels, fitName)
  if (any(vapply(alternatives, function(x) !is.null(x$data), logical(1)))) {
    cli::cli_inform(c(
      "i" = "Alternative-model {.arg data} overrides are recorded but ignored in the current SSE phase."
    ))
  }

  samples <- as.integer(samples)
  projected_total <- samples *
    (length(alternatives) + if (isTRUE(control$estimateSimulation)) 1L else 0L)
  worker_desc <- .workerDescription(control$workers)

  run_dir <- nlmixr2utils::resolveRunDir(
    "sse",
    fitName,
    restart,
    outputDir = outputDir
  )
  if (identical(run_dir$mode, "overwrite") && dir.exists(run_dir$path)) {
    unlink(run_dir$path, recursive = TRUE, force = TRUE)
  }
  dir.create(run_dir$path, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(run_dir$path, "simulations"), showWarnings = FALSE)
  dir.create(file.path(run_dir$path, "fits"), showWarnings = FALSE)

  master_seed <- nlmixr2utils::withRunSeed(
    run_dir$path,
    seed = seed,
    prefix = "sse"
  )
  abs_output_dir <- normalizePath(run_dir$path, mustWork = FALSE)

  cli::cli_rule(left = "SSE")
  cli::cli_inform(c(
    "i" = "Model             : {fitName}",
    "i" = "Parameter source  : {control$parameterSource}",
    "i" = "Replicates        : {samples}",
    "i" = "Projected fits    : {projected_total}",
    "i" = "Execution         : {worker_desc}",
    "i" = "Run directory     : {abs_output_dir}"
  ))
  if (identical(control$parameterSource, "covariance")) {
    cli::cli_inform(c(
      "i" = "Thetas are drawn from {.code fit$cov}; OMEGA and SIGMA remain at point estimates."
    ))
  }
  if (projected_total > 500L) {
    cli::cli_warn(c(
      "!" = "This SSE run projects {projected_total} total fits.",
      "i" = "Adjust {.arg samples} or {.arg alternativeModels} if you want a smaller run."
    ))
  }
  cli::cli_inform(c(
    "i" = "Phase 3 initializes run metadata and output structure; simulation and refit execution arrives in Phase 4."
  ))
  cli::cli_rule()

  timestamp <- Sys.time()
  run_info <- list(
    status = "initialized",
    fitName = fitName,
    samples = samples,
    seed = master_seed,
    parameterSource = control$parameterSource,
    estimateSimulation = control$estimateSimulation,
    projectedTotalFits = projected_total,
    workers = control$workers,
    runDirMode = run_dir$mode,
    alternativeLabels = vapply(alternatives, `[[`, character(1), "label"),
    control = control,
    call = match.call()
  )
  saveRDS(run_info, file.path(abs_output_dir, "run_info.rds"))

  .newNlmixr2SSE(
    runInfo = run_info,
    rawResults = .emptyRawResults(fit, fitName),
    alternativeSpecs = .alternativeSnapshots(alternatives),
    outputDir = abs_output_dir,
    timestamp = timestamp
  )
}
