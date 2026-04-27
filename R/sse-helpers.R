.abortSSE <- function(...) {
  cli::cli_abort(c("!" = ...))
}

.isScalarCharacter <- function(x) {
  is.character(x) && length(x) == 1L && !is.na(x)
}

.validateSSEFit <- function(fit, arg = "fit") {
  if (!inherits(fit, "nlmixr2FitCore")) {
    .abortSSE(paste0("`", arg, "` must inherit from `nlmixr2FitCore`."))
  }

  required <- c("finalUiEnv", "est", "control")
  missing_fields <- required[vapply(
    required,
    function(field) is.null(fit[[field]]),
    logical(1)
  )]
  if (length(missing_fields) > 0L) {
    .abortSSE(
      paste0(
        "`",
        arg,
        "` is missing required field(s): ",
        paste(missing_fields, collapse = ", "),
        " for the current SSE phase."
      )
    )
  }

  invisible(fit)
}

.validateModelLabel <- function(label, arg = "label") {
  if (!.isScalarCharacter(label) || !nzchar(label)) {
    .abortSSE(paste0("`", arg, "` must be a non-empty single string."))
  }
  if (!grepl("^[A-Za-z0-9_.]+$", label)) {
    .abortSSE(
      paste0(
        "`",
        arg,
        "` must contain only letters, numbers, `_`, and `.`."
      )
    )
  }
  invisible(label)
}

.validateFilterInput <- function(filter, arg) {
  tryCatch(
    {
      nlmixr2utils::setupRawResultsFilter(filter)
      invisible(filter)
    },
    error = function(e) {
      .abortSSE(
        paste0(
          "`",
          arg,
          "` is not a valid raw-results filter: ",
          conditionMessage(e)
        )
      )
    }
  )
}

.workerDescription <- function(workers) {
  if (is.null(workers)) {
    return("sequential (using current future::plan())")
  }
  if (identical(workers, "auto")) {
    if (requireNamespace("future", quietly = TRUE)) {
      return(paste0(
        "parallel, auto (",
        future::availableCores(omit = 1L),
        " workers)"
      ))
    }
    return("sequential (future not available)")
  }
  if (identical(as.integer(workers), 1L)) {
    return("sequential (workers = 1)")
  }
  paste0("parallel (", workers, " workers)")
}

.emptyReferenceValues <- function() {
  data.frame(
    parameter = character(0),
    replicate = integer(0),
    value = numeric(0),
    stringsAsFactors = FALSE
  )
}

.emptyInitialValues <- function() {
  data.frame(
    replicate = integer(0),
    parameter = character(0),
    value = numeric(0),
    stringsAsFactors = FALSE
  )
}

.emptyParameterSummary <- function() {
  data.frame(
    model_label = character(0),
    role = character(0),
    parameter = character(0),
    statistic = character(0),
    value = numeric(0),
    n_effective = integer(0),
    matched = logical(0),
    stringsAsFactors = FALSE
  )
}

.emptyOfvSummary <- function() {
  data.frame(
    model_label = character(0),
    statistic = character(0),
    threshold = numeric(0),
    direction = character(0),
    value = numeric(0),
    stringsAsFactors = FALSE
  )
}

.emptyRawResults <- function(fit, model_label) {
  prototype <- nlmixr2utils::rawResultsRow(
    fit = unclass(fit),
    source = "sse",
    hypothesis = "simulation",
    sample = 0L,
    modelLabel = model_label,
    role = "simulation"
  )
  prototype[0, , drop = FALSE]
}

.normalizeAlternativeModels <- function(alternativeModels, simulationLabel) {
  if (is.null(alternativeModels)) {
    return(list())
  }
  if (inherits(alternativeModels, "nlmixr2SSEModel")) {
    alternativeModels <- list(alternativeModels)
  }
  if (
    !is.list(alternativeModels) ||
      any(vapply(
        alternativeModels,
        function(x) !inherits(x, "nlmixr2SSEModel"),
        logical(1)
      ))
  ) {
    .abortSSE(
      "{.arg alternativeModels} must be NULL, an {.cls nlmixr2SSEModel}, or a list of {.cls nlmixr2SSEModel} objects."
    )
  }

  taken <- simulationLabel
  next_auto <- 1L

  for (i in seq_along(alternativeModels)) {
    label <- alternativeModels[[i]]$label
    if (is.null(label)) {
      repeat {
        candidate <- paste0("alt", next_auto)
        next_auto <- next_auto + 1L
        if (!candidate %in% taken) {
          label <- candidate
          break
        }
      }
      alternativeModels[[i]]$label <- label
    }
    if (label %in% taken) {
      .abortSSE(
        paste0(
          "Alternative-model label `",
          label,
          "` collides with the simulation model or another alternative."
        )
      )
    }
    taken <- c(taken, label)
  }

  alternativeModels
}

.alternativeSnapshots <- function(alternativeModels) {
  lapply(alternativeModels, function(spec) {
    list(
      label = spec$label,
      est = spec$est,
      control = spec$control,
      isFit = spec$isFit,
      hasDataOverride = !is.null(spec$data)
    )
  })
}

.validateCovarianceFit <- function(fit) {
  cov_mat <- fit[["cov"]]
  if (
    !is.matrix(cov_mat) ||
      nrow(cov_mat) == 0L ||
      nrow(cov_mat) != ncol(cov_mat)
  ) {
    .abortSSE(
      "{.arg control$parameterSource = \"covariance\"} requires {.arg fit$cov} to be a non-empty square matrix."
    )
  }
  ok <- tryCatch(
    {
      chol(cov_mat)
      TRUE
    },
    error = function(e) FALSE
  )
  if (!ok) {
    .abortSSE(
      "{.arg control$parameterSource = \"covariance\"} requires {.arg fit$cov} to be positive-definite."
    )
  }
  invisible(fit)
}

.newNlmixr2SSE <- function(
  runInfo,
  rawResults,
  alternativeSpecs,
  outputDir,
  timestamp
) {
  structure(
    list(
      runInfo = runInfo,
      referenceValues = .emptyReferenceValues(),
      initialValues = .emptyInitialValues(),
      rawResults = rawResults,
      parameterSummary = .emptyParameterSummary(),
      ofvSummary = .emptyOfvSummary(),
      powerSummary = .emptyOfvSummary(),
      alternativeSpecs = alternativeSpecs,
      outputDir = outputDir,
      timestamp = timestamp
    ),
    class = c("nlmixr2SSE", "list")
  )
}
