#' Configure an `nlmixr2sse` run
#'
#' `runSSEControl()` constructs and validates the control object used by
#' [runSSE()]. Phase 3 implements the argument surface, combination checks, and
#' staged-availability guards so later phases can add the simulation core
#' without changing the external API.
#'
#' @param estimateSimulation Logical. Should the simulation model also be
#'   re-estimated on each simulated dataset?
#' @param refOfv Optional reference OFV override.
#' @param parameterSource Source of simulation parameters: fixed fit values,
#'   canonical raw-results input, or multivariate-normal covariance draws.
#' @param rawresInput Optional raw-results input path or object for
#'   `parameterSource = "rawres"`.
#' @param offsetRawres Integer sample offset for raw-results input.
#' @param inFilter,outFilter Optional raw-results filters. These accept the same
#'   forms as [nlmixr2utils::setupRawResultsFilter()].
#' @param randomEstimationInits Logical. Reserved for raw-results-driven runs.
#' @param updateFix Logical. Reserved for raw-results-driven runs.
#' @param appendColumns Optional character vector of simulated-data columns to
#'   append to SSE outputs.
#' @param simulationPostProcess Reserved for Phase 7. Must be `NULL` in Phase 3.
#' @param initialEtas Reserved for Phase 7. Must be `FALSE` in Phase 3.
#' @param workers Worker setting passed through to shared worker helpers.
#' @param addModels Logical. Add alternative models onto a prior run.
#' @param saveFits,saveDatasets Logical retention flags for later SSE phases.
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
  checkmate::assertFlag(initialEtas)
  checkmate::assertFlag(addModels)
  checkmate::assertFlag(saveFits)
  checkmate::assertFlag(saveDatasets)
  checkmate::assertFlag(overwrite)
  nlmixr2utils::.validateWorkers(workers) # nolint: object_usage_linter.

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
  if (!is.null(simulationPostProcess)) {
    .abortSSE(
      "{.arg simulationPostProcess} is reserved for Phase 7 and must be NULL in the current release."
    )
  }
  if (isTRUE(initialEtas)) {
    .abortSSE(
      "{.arg initialEtas} is reserved for Phase 7 and must remain FALSE in the current release."
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
      addModels = addModels,
      saveFits = saveFits,
      saveDatasets = saveDatasets,
      overwrite = overwrite
    ),
    class = "nlmixr2SSEControl"
  )
}
