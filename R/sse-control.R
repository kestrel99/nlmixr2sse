#' Configure an `nlmixr2sse` run
#'
#' `runSSEControl()` constructs and validates the control object used by
#' [runSSE()]. Phase 7 keeps the Phase 6 uncertainty-mode controls, resumable
#' add-models workflow, raw-results filtering, and estimation-initialization
#' switches, and now also exposes simulation post-processing hooks while still
#' holding `initialEtas` back pending stability work.
#'
#' @param estimateSimulation Logical. Should the simulation model also be
#'   re-estimated on each simulated dataset?
#' @param refOfv Optional reference OFV override.
#' @param parameterSource Source of simulation parameters: fixed fit values,
#'   canonical raw-results input, or multivariate-normal covariance draws.
#' @param rawresInput Optional raw-results input path or object for
#'   `parameterSource = "rawres"`. This is required when
#'   `parameterSource = "rawres"`.
#' @param offsetRawres Integer sample offset for raw-results input.
#' @param inFilter,outFilter Optional raw-results filters. These accept the same
#'   forms as [nlmixr2utils::setupRawResultsFilter()].
#' @param randomEstimationInits Logical. When `TRUE`, raw-results-driven SSE
#'   also uses each selected raw-results row as the starting point for the
#'   estimation step.
#' @param updateFix Logical. When `TRUE`, raw-results-driven SSE updates fixed
#'   model values from the selected raw-results row before refitting.
#' @param appendColumns Optional character vector of simulated-data columns to
#'   append to SSE outputs.
#' @param simulationPostProcess Optional function applied to each simulated
#'   estimation data set after simulated observations have been merged back
#'   onto the original design data. The function must return a data frame. It
#'   may declare any subset of the named arguments `data`, `sample`,
#'   `paramSet`, `solved`, `referenceData`, and `outputDir`.
#' @param initialEtas Reserved for later work. Must remain `FALSE` in the
#'   current release.
#' @param workers Worker setting passed through to shared worker helpers.
#' @param rxThreads Number of rxode2 threads each worker may use. `"auto"`
#'   (the default) divides the machine's cores among the workers. A positive
#'   integer sets the count explicitly; `NULL` defers to rxode2's own default.
#'   Note that rxode2's thread count changes simulated values even under a
#'   fixed seed, so this setting affects reproducibility, not only speed.
#' @param addModels Logical. Extend a completed run by fitting only new
#'   alternative models on the saved simulated datasets.
#' @param saveFits,saveDatasets Logical retention flags used by restart and
#'   add-models workflows.
#' @param overwrite Logical. Reserved output-overwrite flag for later phases.
#'
#' @return An object of class `nlmixr2SSEControl`.
#' @export
runSSEControl <- function(
  estimateSimulation = TRUE,
  refOfv = NULL,
  parameterSource = c("fixed", "rawres", "covariance"),
  rawresInput = NULL,
  offsetRawres = 1L,
  inFilter = NULL,
  outFilter = NULL,
  randomEstimationInits = FALSE,
  updateFix = FALSE,
  appendColumns = NULL,
  simulationPostProcess = NULL,
  initialEtas = FALSE,
  workers = NULL,
  rxThreads = "auto",
  addModels = FALSE,
  saveFits = TRUE,
  saveDatasets = TRUE,
  overwrite = FALSE
) {
  parameterSource <- match.arg(parameterSource)

  checkmate::assertFlag(estimateSimulation)
  if (!is.null(refOfv)) {
    checkmate::assertNumber(refOfv, finite = TRUE)
  }
  checkmate::assertIntegerish(
    offsetRawres,
    len = 1L,
    any.missing = FALSE,
    lower = 0
  )
  checkmate::assertFlag(randomEstimationInits)
  checkmate::assertFlag(updateFix)
  if (!is.null(appendColumns)) {
    checkmate::assertCharacter(
      appendColumns,
      any.missing = FALSE,
      unique = TRUE
    )
  }
  if (!is.null(simulationPostProcess)) {
    checkmate::assertFunction(simulationPostProcess)
  }
  checkmate::assertFlag(initialEtas)
  checkmate::assertFlag(addModels)
  checkmate::assertFlag(saveFits)
  checkmate::assertFlag(saveDatasets)
  checkmate::assertFlag(overwrite)
  nlmixr2utils::.validateWorkers(workers) # nolint: object_usage_linter.
  nlmixr2utils::.validateRxThreads(rxThreads) # nolint: object_usage_linter.

  if (!is.null(inFilter)) {
    .validateFilterInput(inFilter, "inFilter")
  }
  if (!is.null(outFilter)) {
    .validateFilterInput(outFilter, "outFilter")
  }

  if (!is.null(refOfv) && isTRUE(estimateSimulation)) {
    .abortSSE(
      "{.arg refOfv} cannot be supplied when {.arg estimateSimulation} is TRUE."
    )
  }
  if (isTRUE(randomEstimationInits) && parameterSource != "rawres") {
    .abortSSE(
      "{.arg randomEstimationInits} requires {.arg parameterSource = \"rawres\"}."
    )
  }
  if (isTRUE(updateFix) && parameterSource != "rawres") {
    .abortSSE(
      "{.arg updateFix} requires {.arg parameterSource = \"rawres\"}."
    )
  }
  if (!is.null(inFilter) && parameterSource != "rawres") {
    .abortSSE(
      "{.arg inFilter} requires {.arg parameterSource = \"rawres\"}."
    )
  }
  if (!is.null(appendColumns) && isTRUE(addModels)) {
    .abortSSE(
      "{.arg appendColumns} cannot be used together with {.arg addModels}."
    )
  }
  if (isTRUE(initialEtas)) {
    .abortSSE(
      "{.arg initialEtas} is still reserved and must remain FALSE in the current release."
    )
  }

  structure(
    list(
      estimateSimulation = estimateSimulation,
      refOfv = refOfv,
      parameterSource = parameterSource,
      rawresInput = rawresInput,
      offsetRawres = as.integer(offsetRawres),
      inFilter = inFilter,
      outFilter = outFilter,
      randomEstimationInits = randomEstimationInits,
      updateFix = updateFix,
      appendColumns = appendColumns,
      simulationPostProcess = simulationPostProcess,
      initialEtas = initialEtas,
      workers = workers,
      rxThreads = rxThreads,
      addModels = addModels,
      saveFits = saveFits,
      saveDatasets = saveDatasets,
      overwrite = overwrite
    ),
    class = "nlmixr2SSEControl"
  )
}
