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
#' @param covarianceDraw How `parameterSource = "covariance"` draws parameters.
#'   `"independent_iw"` (the default) draws THETA multivariate-Normal and OMEGA
#'   from a mean-centred inverse-Wishart per OMEGA block, independently of each
#'   other. `"joint"` instead draws THETA and OMEGA together from `fit$cov` on a
#'   log-Cholesky-transformed scale, incorporating the estimated THETA/OMEGA
#'   covariance through a first-order delta approximation. Both give
#'   positive-definite OMEGA draws.
#'
#'   `"joint"` uses information `"independent_iw"` discards, but its non-linear
#'   back-transform inflates OMEGA means by roughly `exp(SE^2 / (2 * Omega^2))`
#'   — about 2% at 20% relative standard error, but 65% at 100% — so it is
#'   opt-in rather than the default. Neither mode implements NONMEM's
#'   `$PRIOR NWPRI`, and neither claims PsN parity.
#' @param omegaRseWarn Relative-standard-error threshold above which a drawn
#'   OMEGA variance triggers a weak-identification warning naming the block and
#'   eta. Defaults to `0.5`.
#' @param rawresInput Optional raw-results input path or object for
#'   `parameterSource = "rawres"`. This is required when
#'   `parameterSource = "rawres"`.
#' @param offsetRawres Integer sample offset for raw-results input.
#' @param inFilter,outFilter Optional raw-results filters. These accept the same
#'   forms as [nlmixr2utils::setupRawResultsFilter()].
#' @param referenceInitials,alternativeInitials Estimation starting-value
#'   policy for, respectively, the simulation (reference) model refit and the
#'   alternative-model refits. `"model"` (the default for both) starts each
#'   fit from the model's own initial estimates, as before. `"simulation"`
#'   starts a fit from that replicate's generating parameter vector instead --
#'   a numerical intervention only: it changes optimiser starting values, not
#'   the simulated data. It only has an effect when
#'   `parameterSource = "rawres"` -- the same restriction the deprecated
#'   `randomEstimationInits` enforced as a hard error; setting either argument
#'   to `"simulation"` under another mode instead warns (not errors, since
#'   setting one role while the other keeps its default is legitimate). The
#'   two other modes are inert for different reasons: under `"fixed"` every
#'   replicate shares one generating vector, so simulation and model starts
#'   are byte-identical -- a genuine no-op. Under `"covariance"` each
#'   replicate draws a genuinely different vector, so a simulation start
#'   *would* differ from a model start, but `.fitTaskRecord()` only ever
#'   applies "simulation" starts when `parameterSource` is `"rawres"` -- this
#'   is unimplemented, not inert. Splitting the setting by role lets a
#'   reference-only or alternative-only sensitivity study be run without
#'   touching the other role's starts.
#' @param randomEstimationInits Deprecated. Use `referenceInitials` and
#'   `alternativeInitials` instead. `TRUE` maps to
#'   `referenceInitials = "simulation"` and `alternativeInitials =
#'   "simulation"`; `FALSE` maps both to `"model"`.
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
  covarianceDraw = c("independent_iw", "joint"),
  omegaRseWarn = 0.5,
  rawresInput = NULL,
  offsetRawres = 1L,
  inFilter = NULL,
  outFilter = NULL,
  referenceInitials = c("model", "simulation"),
  alternativeInitials = c("model", "simulation"),
  randomEstimationInits = lifecycle::deprecated(),
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
  # Captured before match.arg() reassigns the arguments below -- missing()
  # inspects how the argument was bound at call time, and reassigning the
  # variable replaces that binding, so a missing() call after match.arg()
  # would always report FALSE regardless of what the caller supplied.
  referenceInitialsMissing <- missing(referenceInitials)
  alternativeInitialsMissing <- missing(alternativeInitials)

  parameterSource <- match.arg(parameterSource)
  covarianceDrawMissing <- missing(covarianceDraw)
  covarianceDraw <- match.arg(covarianceDraw)
  omegaRseWarnMissing <- missing(omegaRseWarn)
  # match.arg()'s own error reads "should be one of", which does not match
  # this package's "{.arg x} must be one of ..." phrasing used elsewhere for
  # invalid enums, so its error is caught and rephrased. Explicit `choices`
  # keeps match.arg's partial-matching behavior working even though the call
  # is no longer a direct, unwrapped reference to the formal argument.
  referenceInitials <- tryCatch(
    match.arg(referenceInitials, c("model", "simulation")),
    error = function(e) {
      .abortSSE(
        "{.arg referenceInitials} must be one of {.val {c('model', 'simulation')}}."
      )
    }
  )
  alternativeInitials <- tryCatch(
    match.arg(alternativeInitials, c("model", "simulation")),
    error = function(e) {
      .abortSSE(
        "{.arg alternativeInitials} must be one of {.val {c('model', 'simulation')}}."
      )
    }
  )

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
  checkmate::assertFlag(updateFix)
  checkmate::assertNumber(
    omegaRseWarn,
    lower = 0,
    finite = TRUE,
    .var.name = "omegaRseWarn"
  )
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

  # randomEstimationInits applied its single TRUE/FALSE value to every model
  # role at once, and only ever meant anything under parameterSource =
  # "rawres" (every other mode gives every replicate the same generating
  # vector, so "start at the generating value" is just "start at the fitted
  # value" for every fit) -- so it hard-errored outside that mode. Preserve
  # that exact behavior for the deprecated argument. referenceInitials /
  # alternativeInitials generalize the same "simulation" policy to modes
  # where it is merely inert rather than nonsensical (a fixed or covariance
  # run still resolves a real, if uniform-per-replicate, generating vector),
  # so setting them outside rawres mode is left as a silent no-op --
  # .fitTaskRecord() only ever applies "simulation" starts when
  # parameterSource is "rawres" -- rather than an error.
  if (lifecycle::is_present(randomEstimationInits)) {
    lifecycle::deprecate_soft(
      "0.1",
      "runSSEControl(randomEstimationInits)",
      details = "Use `referenceInitials` and `alternativeInitials` instead."
    )
    checkmate::assertFlag(randomEstimationInits)
    mapped <- if (isTRUE(randomEstimationInits)) "simulation" else "model"
    # Contradiction is checked ahead of the rawres restriction below: a call
    # that both contradicts itself and violates that restriction should
    # report the contradiction, since it is the more fundamental problem.
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
    if (isTRUE(randomEstimationInits) && parameterSource != "rawres") {
      .abortSSE(
        "{.arg randomEstimationInits} requires {.arg parameterSource = \"rawres\"}."
      )
    }
    referenceInitials <- mapped
    alternativeInitials <- mapped
  }

  # Unlike the deprecated flag, an explicit "simulation" outside "rawres"
  # mode warns rather than errors -- setting one role while the other keeps
  # its default "model" is legitimate, so this cannot be a hard error. But
  # every other mode-gated argument in this function (updateFix, inFilter,
  # covarianceDraw, omegaRseWarn) never silently swallows an
  # explicitly-set-but-inapplicable option, so this cannot be silent either.
  if (identical(referenceInitials, "simulation") && parameterSource != "rawres") {
    cli::cli_warn(c(
      "!" = "{.arg referenceInitials = \"simulation\"} has no effect unless {.arg parameterSource = \"rawres\"}.",
      "i" = "Under {.val fixed} every replicate shares one generating vector, so the two policies coincide.",
      "i" = "Under {.val covariance} the vectors differ per replicate, but simulation starts are not yet wired up."
    ))
  }
  if (identical(alternativeInitials, "simulation") && parameterSource != "rawres") {
    cli::cli_warn(c(
      "!" = "{.arg alternativeInitials = \"simulation\"} has no effect unless {.arg parameterSource = \"rawres\"}.",
      "i" = "Under {.val fixed} every replicate shares one generating vector, so the two policies coincide.",
      "i" = "Under {.val covariance} the vectors differ per replicate, but simulation starts are not yet wired up."
    ))
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
  if (!covarianceDrawMissing && parameterSource != "covariance") {
    .abortSSE(
      "{.arg covarianceDraw} requires {.arg parameterSource = \"covariance\"}."
    )
  }
  if (!omegaRseWarnMissing && parameterSource != "covariance") {
    .abortSSE(
      "{.arg omegaRseWarn} requires {.arg parameterSource = \"covariance\"}."
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
      covarianceDraw = covarianceDraw,
      omegaRseWarn = omegaRseWarn,
      rawresInput = rawresInput,
      offsetRawres = as.integer(offsetRawres),
      inFilter = inFilter,
      outFilter = outFilter,
      referenceInitials = referenceInitials,
      alternativeInitials = alternativeInitials,
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
