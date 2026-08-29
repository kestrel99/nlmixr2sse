.abortSSE <- function(..., .envir = parent.frame()) {
  cli::cli_abort(c("!" = ...), .envir = .envir)
}

.isScalarCharacter <- function(x) {
  is.character(x) && length(x) == 1L && !is.na(x)
}

.truncateSSEMessage <- function(message, width = 500L) {
  message <- paste(as.character(message), collapse = " ")
  message <- gsub("[\r\n\t]+", " ", message)
  message <- gsub(" +", " ", message)
  message <- trimws(message)
  if (nchar(message, type = "chars") <= width) {
    return(message)
  }
  paste0(substr(message, 1L, width - 3L), "...")
}

.fitField <- function(fit, field) {
  value <- tryCatch(fit[[field]], error = function(e) NULL)
  if (!is.null(value)) {
    return(value)
  }

  obj <- fit
  tryCatch(
    eval(call("$", quote(obj), as.name(field))),
    error = function(e) NULL
  )
}

.validateSSEFit <- function(fit, arg = "fit") {
  if (!inherits(fit, "nlmixr2FitCore")) {
    .abortSSE(paste0("`", arg, "` must inherit from `nlmixr2FitCore`."))
  }

  required <- c("est", "control")
  missing_fields <- required[vapply(
    required,
    function(field) is.null(.fitField(fit, field)),
    logical(1)
  )]
  has_ui <- inherits(.fitField(fit, "ui"), "rxUi") ||
    inherits(.fitField(fit, "finalUiEnv"), "rxUi")
  if (!has_ui) {
    missing_fields <- c(missing_fields, "ui/finalUiEnv")
  }
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

.workerDescription <- function(workers, rxThreads = NULL) {
  base <- if (is.null(workers)) {
    "sequential (using current future::plan())"
  } else if (identical(workers, "auto")) {
    if (requireNamespace("future", quietly = TRUE)) {
      paste0(
        "parallel, auto (",
        future::availableCores(omit = 1L),
        " workers)"
      )
    } else {
      "sequential (future not available)"
    }
  } else if (identical(as.integer(workers), 1L)) {
    "sequential (workers = 1)"
  } else {
    paste0("parallel (", workers, " workers)")
  }

  if (is.null(rxThreads)) {
    return(base)
  }
  paste0(
    base,
    ", ",
    rxThreads,
    " rxode2 thread",
    if (identical(as.integer(rxThreads), 1L)) "" else "s",
    "/worker"
  )
}

#' Choose the rxThreads value to record in run_info
#'
#' An add-models run refits new alternatives against the SAVED simulated
#' datasets, so the recorded thread count should keep describing how those
#' datasets were produced rather than being overwritten by whatever the
#' add-models call happens to resolve to.
#' @noRd
.addModelsRxThreads <- function(
  existingRunInfo,
  resolvedRxThreads,
  addModels
) {
  if (!isTRUE(addModels)) {
    return(resolvedRxThreads)
  }
  recorded <- existingRunInfo$rxThreads
  if (is.null(recorded)) {
    return(resolvedRxThreads)
  }
  as.integer(recorded)
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
    fit = fit,
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

.schemaSnapshot <- function(schema) {
  list(
    thetaCols = unname(schema$thetaCols %||% character(0)),
    omegaCols = unname(schema$omegaCols %||% character(0)),
    sigmaCols = unname(schema$sigmaCols %||% character(0))
  )
}

.fitSpecsSnapshot <- function(fitSpecs) {
  lapply(fitSpecs, function(spec) {
    list(
      label = spec$label,
      role = spec$role,
      hypothesis = spec$hypothesis,
      est = spec$est,
      schema = .schemaSnapshot(spec$schema)
    )
  })
}

.rawResultsSchemaFromHeader <- function(header) {
  list(
    columns = unname(header$columns %||% character(0)),
    baseCols = unname(header$base_cols %||% character(0)),
    thetaCols = unname(header$theta_cols %||% character(0)),
    omegaCols = unname(header$omega_cols %||% character(0)),
    sigmaCols = unname(header$sigma_cols %||% character(0)),
    seCols = unname(header$se_cols %||% character(0)),
    schemaVersion = as.integer(header$schema_version %||% 1L)
  )
}

.fitSpecsFromRunInfo <- function(runInfo, rawResults = NULL) {
  if (is.list(runInfo$fitSpecs) && length(runInfo$fitSpecs) > 0L) {
    return(runInfo$fitSpecs)
  }

  if (is.null(rawResults) || nrow(rawResults) == 0L) {
    return(list())
  }

  header <- attr(rawResults, "rawResultsHeader", exact = TRUE)
  if (is.null(header)) {
    .abortSSE(
      "Completed raw results are missing the canonical header metadata."
    )
  }
  schema <- .rawResultsSchemaFromHeader(header)
  labels <- unique(rawResults$model_label)
  lapply(labels, function(label) {
    rows <- rawResults$model_label == label
    list(
      label = label,
      role = rawResults$role[match(TRUE, rows)],
      hypothesis = rawResults$hypothesis[match(TRUE, rows)],
      est = NA_character_,
      schema = .schemaSnapshot(schema)
    )
  })
}

.specMapFromSnapshots <- function(fitSpecsSnapshot) {
  if (length(fitSpecsSnapshot) == 0L) {
    return(list())
  }
  stats::setNames(
    lapply(fitSpecsSnapshot, function(spec) {
      spec$schema
    }),
    vapply(fitSpecsSnapshot, `[[`, character(1), "label")
  )
}

.alternativeSnapshotsFromRunInfo <- function(runInfo) {
  if (is.list(runInfo$alternativeSpecs)) {
    return(runInfo$alternativeSpecs)
  }

  fit_specs <- .fitSpecsFromRunInfo(runInfo)
  lapply(
    Filter(function(spec) identical(spec$role, "alternative"), fit_specs),
    function(spec) {
      list(
        label = spec$label,
        est = spec$est,
        control = NULL,
        isFit = NA,
        hasDataOverride = FALSE
      )
    }
  )
}

.validateCovarianceFit <- function(fit) {
  cov_mat <- .fitField(fit, "cov")
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

.thetaVectorFromFit <- function(fit) {
  theta <- .fitField(fit, "theta")
  if (is.null(theta) || length(theta) == 0L) {
    return(stats::setNames(numeric(0), character(0)))
  }
  theta_vals <- as.numeric(theta)
  names(theta_vals) <- names(theta)
  theta_vals
}

.paramSetFromFit <- function(fit) {
  list(
    theta = .thetaVectorFromFit(fit),
    omega = if (is.matrix(.fitField(fit, "omega"))) {
      .fitField(fit, "omega")
    } else {
      matrix(numeric(0), 0L, 0L)
    },
    sigma = if (is.matrix(.fitField(fit, "sigma"))) {
      .fitField(fit, "sigma")
    } else {
      matrix(numeric(0), 0L, 0L)
    }
  )
}

.getFitUi <- function(fit) {
  ui <- .fitField(fit, "ui")
  if (inherits(ui, "rxUi")) {
    return(ui)
  }
  ui <- .fitField(fit, "finalUiEnv")
  if (inherits(ui, "rxUi")) {
    return(ui)
  }
  .abortSSE("The supplied fit does not expose a usable rxUi simulation model.")
}

.getFitData <- function(fit) {
  data <- tryCatch(
    nlme::getData(fit),
    error = function(e) NULL
  )
  if (is.null(data)) {
    data <- fit[["data"]]
  }
  if (is.null(data)) {
    .abortSSE(
      "Could not recover the original estimation data from the reference fit."
    )
  }
  as.data.frame(data, stringsAsFactors = FALSE)
}

.inferStudySampleInfo <- function(data) {
  checkmate::assertDataFrame(data, min.rows = 1L)

  id_match <- match("id", tolower(names(data)))
  if (!is.na(id_match)) {
    id_col <- names(data)[[id_match]]
    ids <- data[[id_col]]
    ids <- ids[!is.na(ids)]
    return(list(
      studySampleSize = as.integer(length(unique(ids))),
      studySampleUnit = "subjects",
      studyIdColumn = id_col,
      studyObservationCount = as.integer(nrow(data))
    ))
  }

  list(
    studySampleSize = as.integer(nrow(data)),
    studySampleUnit = "observations",
    studyIdColumn = NA_character_,
    studyObservationCount = as.integer(nrow(data))
  )
}

.rawResultsSchemaForFit <- function(fit) {
  # Access each field through the fit's own dispatch before building the plain
  # list that rawResultsSchema can consume safely. unclass() strips nlmixr2's
  # active bindings and method dispatch, causing theta/omega/sigma to return
  # NULL or unnamed values and producing a schema with no parameter columns.
  fit_like <- list(
    theta = tryCatch(fit$theta, error = function(e) NULL),
    omega = tryCatch(fit$omega, error = function(e) NULL),
    sigma = tryCatch(fit$sigma, error = function(e) NULL),
    iniDf = tryCatch(fit$iniDf, error = function(e) NULL),
    parFixedDf = tryCatch(fit$parFixedDf, error = function(e) NULL)
  )
  nlmixr2utils::rawResultsSchema(fit_like)
}

.schemaParamCols <- function(schema) {
  c(schema$thetaCols, schema$omegaCols, schema$sigmaCols)
}

.schemaParamPresence <- function(schema) {
  cols <- .schemaParamCols(schema)
  stats::setNames(rep(TRUE, length(cols)), cols)
}

.fitLikeFromUi <- function(ui) {
  ini_df <- ui[["iniDf"]]
  if (!is.data.frame(ini_df)) {
    return(list(
      theta = numeric(0),
      omega = matrix(numeric(0), 0L, 0L),
      sigma = matrix(numeric(0), 0L, 0L),
      iniDf = data.frame(stringsAsFactors = FALSE)
    ))
  }

  theta_rows <- !is.na(ini_df$ntheta)
  if ("fix" %in% names(ini_df)) {
    theta_rows <- theta_rows & !ini_df$fix
  }
  theta_names <- ini_df$name[theta_rows]
  theta_vals <- if ("est" %in% names(ini_df)) {
    as.numeric(ini_df$est[theta_rows])
  } else {
    rep(NA_real_, sum(theta_rows))
  }
  keep_theta <- !is.na(theta_names) & nzchar(theta_names)
  theta <- stats::setNames(theta_vals[keep_theta], theta_names[keep_theta])

  omega_rows <- ini_df[!is.na(ini_df$neta1), , drop = FALSE]
  if ("fix" %in% names(omega_rows)) {
    omega_rows <- omega_rows[!omega_rows$fix, , drop = FALSE]
  }

  if (nrow(omega_rows) == 0L) {
    omega <- matrix(numeric(0), 0L, 0L)
  } else {
    n_eta <- max(
      c(0L, as.integer(omega_rows$neta1), as.integer(omega_rows$neta2)),
      na.rm = TRUE
    )
    omega <- matrix(0, n_eta, n_eta)
    diag_rows <- omega_rows[
      omega_rows$neta1 == omega_rows$neta2,
      ,
      drop = FALSE
    ]
    eta_names <- as.character(seq_len(n_eta))
    if (nrow(diag_rows) > 0L) {
      eta_names[diag_rows$neta1] <- diag_rows$name
    }
    dimnames(omega) <- list(eta_names, eta_names)

    est_vals <- if ("est" %in% names(omega_rows)) {
      as.numeric(omega_rows$est)
    } else {
      rep(0, nrow(omega_rows))
    }
    for (i in seq_len(nrow(omega_rows))) {
      r <- as.integer(omega_rows$neta1[[i]])
      c_idx <- as.integer(omega_rows$neta2[[i]])
      omega[r, c_idx] <- est_vals[[i]]
      omega[c_idx, r] <- est_vals[[i]]
    }
  }

  list(
    theta = theta,
    omega = omega,
    sigma = matrix(numeric(0), 0L, 0L),
    iniDf = ini_df
  )
}

.schemaFromAlternativeSpec <- function(spec) {
  if (isTRUE(spec$isFit)) {
    return(.rawResultsSchemaForFit(spec$model))
  }

  ui <- spec$ui
  if (!inherits(ui, "rxUi")) {
    ui <- tryCatch(
      nlmixr2est::nlmixr2(ui),
      error = function(e) {
        .abortSSE(
          "Could not parse alternative model {.val {spec$label}} into an rxUi object: {conditionMessage(e)}"
        )
      }
    )
  }
  nlmixr2utils::rawResultsSchema(.fitLikeFromUi(ui))
}

.mergeRawResultsSchemas <- function(schemas) {
  theta_cols <- unique(unlist(
    lapply(schemas, `[[`, "thetaCols"),
    use.names = FALSE
  ))
  omega_cols <- unique(unlist(
    lapply(schemas, `[[`, "omegaCols"),
    use.names = FALSE
  ))
  sigma_cols <- unique(unlist(
    lapply(schemas, `[[`, "sigmaCols"),
    use.names = FALSE
  ))
  param_cols <- c(theta_cols, omega_cols, sigma_cols)

  list(
    columns = c(
      c(
        "source",
        "hypothesis",
        "sample",
        "model_label",
        "role",
        "minimization_successful",
        "covariance_step_successful",
        "estimate_near_boundary",
        "significant_digits",
        "condition_number",
        "objf",
        "error_message"
      ),
      param_cols,
      paste0(param_cols, ".se")
    ),
    baseCols = c(
      "source",
      "hypothesis",
      "sample",
      "model_label",
      "role",
      "minimization_successful",
      "covariance_step_successful",
      "estimate_near_boundary",
      "significant_digits",
      "condition_number",
      "objf",
      "error_message"
    ),
    thetaCols = theta_cols,
    omegaCols = omega_cols,
    sigmaCols = sigma_cols,
    seCols = paste0(param_cols, ".se"),
    schemaVersion = 1L
  )
}

.ensureRxodeCache <- function() {
  if (!requireNamespace("rxode2", quietly = TRUE)) {
    return(invisible(NULL))
  }
  try(rxode2::rxCreateCache(), silent = TRUE)
  invisible(NULL)
}

.simulationModelSpec <- function(fit, fitName) {
  list(
    label = fitName,
    ui = .getFitUi(fit),
    est = as.character(.fitField(fit, "est")),
    control = .fitField(fit, "control"),
    isFit = TRUE,
    role = "simulation",
    hypothesis = "simulation",
    schema = .rawResultsSchemaForFit(fit)
  )
}

.fitTaskSpecs <- function(fit, fitName, alternatives, control) {
  specs <- list()
  if (isTRUE(control$estimateSimulation)) {
    specs[[1L]] <- .simulationModelSpec(fit, fitName)
  }
  if (length(alternatives) > 0L) {
    alt_specs <- lapply(seq_along(alternatives), function(i) {
      spec <- alternatives[[i]]
      list(
        label = spec$label,
        ui = spec$ui,
        est = spec$est,
        control = spec$control,
        isFit = spec$isFit,
        role = "alternative",
        hypothesis = paste0("alternative_", i),
        schema = .schemaFromAlternativeSpec(spec)
      )
    })
    specs <- c(specs, alt_specs)
  }
  specs
}

.simulationParamRecord <- function(paramSet, schema) {
  row <- nlmixr2utils::rawResultsRow(
    fit = NULL,
    source = "sse",
    hypothesis = "simulation",
    sample = 1L,
    modelLabel = "simulation",
    role = "simulation",
    theta = paramSet$theta,
    omega = paramSet$omega,
    sigma = paramSet$sigma,
    schema = schema
  )
  as.numeric(row[1L, .schemaParamCols(schema), drop = TRUE])
}

.simulationRowIdCol <- ".sse_row_id"

.normalizeRawresInput <- function(rawresInput, inFilter = NULL) {
  rawres <- if (is.character(rawresInput) && length(rawresInput) == 1L) {
    nlmixr2utils::readRawResults(rawresInput)
  } else {
    rawresInput
  }
  if (!is.data.frame(rawres)) {
    .abortSSE(
      "{.arg control$rawresInput} must resolve to a canonical raw-results data frame or path."
    )
  }

  if (
    is.null(inFilter) &&
      "source" %in% names(rawres) &&
      "sample" %in% names(rawres) &&
      anyDuplicated(rawres$sample) > 0L
  ) {
    sim_rows <- rawres$source == "sse" &
      (("role" %in% names(rawres) & rawres$role == "simulation") |
        ("hypothesis" %in% names(rawres) & rawres$hypothesis == "simulation"))
    sim_rawres <- rawres[sim_rows, , drop = FALSE]
    if (nrow(sim_rawres) > 0L && anyDuplicated(sim_rawres$sample) == 0L) {
      rawres <- sim_rawres
    }
  }

  rawres
}

.resolveRawresParameterSets <- function(
  fit,
  samples,
  control
) {
  if (is.null(control$rawresInput)) {
    .abortSSE(
      "{.arg control$rawresInput} is required when {.arg control$parameterSource = \"rawres\"}."
    )
  }

  rawres <- .normalizeRawresInput(
    control$rawresInput,
    inFilter = control$inFilter
  )
  params <- nlmixr2utils::parseRawResultsParams(
    rawres,
    fit,
    offset = control$offsetRawres,
    filter = control$inFilter
  )

  if (length(params) < samples) {
    if (!is.null(control$inFilter)) {
      .abortSSE(
        "Only {length(params)} raw-results parameter set{?s} remain after applying {.arg control$offsetRawres} and {.arg control$inFilter}, but {.arg samples = {samples}} requires more."
      )
    }
    .abortSSE(
      "Only {length(params)} raw-results parameter set{?s} remain after applying {.arg control$offsetRawres}, but {.arg samples = {samples}} requires more."
    )
  }

  params <- params[seq_len(samples)]
  records <- lapply(seq_along(params), function(i) {
    list(
      sample = as.integer(i),
      sourceSample = as.integer(params[[i]]$sample),
      source = params[[i]]$source,
      hypothesis = params[[i]]$hypothesis,
      modelLabel = params[[i]]$modelLabel,
      role = params[[i]]$role,
      theta = params[[i]]$theta,
      omega = params[[i]]$omega,
      sigma = params[[i]]$sigma
    )
  })

  input_path <- attr(rawres, "rawResultsPath", exact = TRUE)
  list(
    records = records,
    info = list(
      mode = "rawres",
      input = if (is.null(input_path)) {
        "<raw-results-data.frame>"
      } else {
        input_path
      },
      offsetRawres = control$offsetRawres,
      inFilter = control$inFilter,
      sourceSamples = vapply(records, `[[`, integer(1), "sourceSample"),
      estimationInitialValues = if (isTRUE(control$randomEstimationInits)) {
        "rawres"
      } else {
        "reference_fit"
      },
      updateFixedValues = isTRUE(control$updateFix)
    )
  )
}

.alignedCovariance <- function(fit) {
  .validateCovarianceFit(fit)
  theta <- .thetaVectorFromFit(fit)
  cov_mat <- .fitField(fit, "cov")
  cov_names <- rownames(cov_mat)

  if (is.null(cov_names) || is.null(colnames(cov_mat))) {
    if (length(theta) != nrow(cov_mat)) {
      .abortSSE(
        "{.arg fit$cov} must have row and column names when it does not cover every theta parameter."
      )
    }
    cov_names <- names(theta)
    dimnames(cov_mat) <- list(cov_names, cov_names)
  }

  if (!identical(rownames(cov_mat), colnames(cov_mat))) {
    .abortSSE(
      "{.arg fit$cov} row and column names must match for covariance-driven SSE."
    )
  }
  if (anyDuplicated(cov_names) > 0L || any(!nzchar(cov_names))) {
    .abortSSE(
      "{.arg fit$cov} must use unique, non-empty theta names."
    )
  }
  # fit$cov now carries OMEGA entries (om.<eta> / cov.<eta>.<eta>) alongside the
  # thetas. Partition them: thetas are drawn multivariate-Normal from their own
  # sub-block, OMEGA is drawn separately from an inverse-Wishart, and the
  # THETA<->OMEGA cross-terms are deliberately discarded (this mode shares NWPRI's independence factorization, though it is not
  # NWPRI -- the OMEGA density differs).
  # .fitField(), not fit[["ui"]]: real nlmixr2 fits expose $ via a custom S3
  # method but have NO [[ method, so fit[["ui"]] falls through to [[.data.frame
  # and silently returns NULL -- yielding an empty OMEGA table and dropping
  # every OMEGA entry. .fitField() tries [[ then falls back to $ dispatch.
  omega_entries <- .omegaEntryTable(.uiOmegaInfo(.fitField(fit, "ui")))
  omega_names <- omega_entries$covName

  theta_names <- intersect(cov_names, names(theta))
  unknown <- setdiff(cov_names, c(theta_names, omega_names))
  if (length(unknown) > 0L) {
    .abortSSE(
      "{.arg fit$cov} contains name{?s} {.val {unknown}} that match neither {.arg fit$theta} nor an OMEGA entry."
    )
  }

  # per-eta standard errors, from the diagonal OMEGA entries present in fit$cov
  omega_mat <- .fitField(fit, "omega")
  eta_names <- if (is.matrix(omega_mat)) rownames(omega_mat) else character(0)
  omega_se <- stats::setNames(
    rep(NA_real_, length(eta_names)),
    eta_names
  )
  diag_entries <- omega_entries[omega_entries$diagonal, , drop = FALSE]
  for (i in seq_len(nrow(diag_entries))) {
    nm <- diag_entries$covName[[i]]
    eta <- diag_entries$rowName[[i]]
    if (isTRUE(diag_entries$fix[[i]])) {
      next
    }
    if (nm %in% cov_names && eta %in% eta_names) {
      omega_se[[eta]] <- sqrt(cov_mat[nm, nm])
    }
  }

  list(
    theta = theta,
    cov = cov_mat[theta_names, theta_names, drop = FALSE],
    # the whole matrix, kept so the "joint" mode can pull a theta+omega
    # sub-block including the cross-terms that "independent_iw" discards
    fullCov = cov_mat,
    drawNames = theta_names,
    omegaSe = omega_se,
    omegaEntries = omega_entries
  )
}

.resolveCovarianceParameterSets <- function(
  fit,
  samples,
  outputDir,
  schema,
  control
) {
  # No fallback default here: runSSEControl() already resolves covarianceDraw
  # via match.arg(), so a NULL means the control was built by something other
  # than runSSEControl() and we should not silently invent a mode.
  draw_mode <- control$covarianceDraw
  aligned <- .alignedCovariance(fit)
  base_params <- .paramSetFromFit(fit)

  has_theta_draw <- length(aligned$drawNames) > 0L
  chol_cov <- if (has_theta_draw) chol(aligned$cov) else NULL

  omega0 <- base_params$omega
  n_eta <- if (is.matrix(omega0)) nrow(omega0) else 0L
  blocks <- .omegaBlocks(aligned$omegaEntries, n_eta)
  omega_se <- aligned$omegaSe

  # which omega entries actually vary: those in a block with at least one
  # usable standard error
  # Coverage policy: a block is drawn ONLY when every declared entry in it is
  # unfixed and present in fit$cov. Both modes consume this same result -- do
  # NOT re-derive drawability from standard errors here.
  coverage <- .drawableOmegaBlocks(
    blocks,
    aligned$omegaEntries,
    rownames(aligned$fullCov)
  )
  drawn_blocks <- coverage$drawable

  # report exactly the entries of the drawn blocks, so the partition matches
  # what actually varies
  entriesOf <- function(idx) {
    inBlock <- aligned$omegaEntries$row %in% idx &
      aligned$omegaEntries$col %in% idx
    .matrixCoordLabelsForEntries(
      aligned$omegaEntries[inBlock, , drop = FALSE]
    )
  }
  drawn_omega_cols <- unlist(
    lapply(drawn_blocks, entriesOf),
    use.names = FALSE
  ) %||% character(0)

  # surface why any block was held, so the run record explains itself
  for (h in coverage$held) {
    cli::cli_inform(c(
      "i" = "OMEGA block {.val {paste(rownames(omega0)[h$index], collapse = ', ')}} held at its fitted values ({h$reason})."
    ))
  }

  # Warn once per run per drawable block about weakly identified variances --
  # NOT once per replicate, which would emit `samples` identical warnings.
  for (idx in drawn_blocks) {
    .warnWeakOmega(
      omega0[idx, idx, drop = FALSE],
      omega_se[idx],
      idx,
      control$omegaRseWarn,
      draw_mode
    )
  }

  # "joint" needs one spec built once: the transformed covariance over
  # c(theta, stacked omega block vectors), restricted to the drawn entries.
  # Build the joint spec whenever there is ANYTHING to draw -- covered thetas,
  # drawable OMEGA blocks, or both. Gating on OMEGA alone would silently stop
  # drawing thetas for a theta-only fit$cov (e.g. foceiControl(covFull =
  # FALSE)), which the coverage policy explicitly supports.
  # .jointDrawSpec() handles an empty `blocks` list: the Jacobian is NULL, so
  # B collapses to diag(nTheta) and Sigma_T is just the theta covariance.
  joint_spec <- NULL
  if (
    identical(draw_mode, "joint") &&
      (has_theta_draw || length(drawn_blocks) > 0L)
  ) {
    joint_names <- c(
      aligned$drawNames,
      unlist(lapply(drawn_blocks, function(idx) {
        inBlock <- aligned$omegaEntries$row %in% idx &
          aligned$omegaEntries$col %in% idx
        aligned$omegaEntries$covName[inBlock]
      }), use.names = FALSE)
    )
    # Defensive only: .drawableOmegaBlocks() already guarantees every entry of
    # a drawable block is present. Reaching this means the coverage policy and
    # this call disagree, which is an internal inconsistency, not user error.
    missing_joint <- setdiff(joint_names, rownames(aligned$fullCov))
    if (length(missing_joint) > 0L) {
      .abortSSE(
        "Internal error: OMEGA entr{?y/ies} {.val {missing_joint}} passed the coverage check but {?is/are} absent from {.arg fit$cov}."
      )
    }
    joint_spec <- .jointDrawSpec(
      theta = aligned$theta[aligned$drawNames],
      blocks = lapply(drawn_blocks, function(idx) {
        list(omega = omega0[idx, idx, drop = FALSE], index = idx)
      }),
      sigma = aligned$fullCov[joint_names, joint_names, drop = FALSE]
    )
  }

  records <- lapply(seq_len(samples), function(sample_id) {
    drawn <- nlmixr2utils::withRunSeed(
      outputDir,
      key = paste0("covariance-", sample_id),
      prefix = "sse",
      expr = {
        if (identical(draw_mode, "joint")) {
          if (is.null(joint_spec)) {
            # nothing is drawable at all: no covered thetas AND no drawable
            # OMEGA block. Everything keeps its fitted value.
            list(theta = numeric(0), omega = omega0)
          } else {
            # THETA and OMEGA drawn TOGETHER, incorporating their covariance.
            j <- .drawJoint(joint_spec)
            om <- omega0
            for (k in seq_along(drawn_blocks)) {
              idx <- drawn_blocks[[k]]
              om[idx, idx] <- j$omega[[k]]
            }
            list(theta = unname(j$theta), omega = om)
          }
        } else {
          # "independent_iw": THETA and OMEGA drawn INDEPENDENTLY. This shares NWPRI's
          # broad independence factorization, but is NOT NWPRI -- the OMEGA
          # density differs (see the mode comparison in the design doc).
          theta_draw <- if (has_theta_draw) {
            noise <- stats::rnorm(length(aligned$drawNames))
            aligned$theta[aligned$drawNames] + as.numeric(noise %*% chol_cov)
          } else {
            numeric(0)
          }
          omega_draw <- if (n_eta > 0L) {
            # drawn_blocks, NOT blocks: .drawOmega() requires the approved
            # drawable list. Passing the raw block list would redraw blocks the
            # coverage policy held fixed, mutating fixed or uncovered entries.
            .drawOmega(omega0, drawn_blocks, omega_se)
          } else {
            omega0
          }
          list(theta = theta_draw, omega = omega_draw)
        }
      }
    )

    theta <- aligned$theta
    if (has_theta_draw && length(drawn$theta) > 0L) {
      theta[aligned$drawNames] <- drawn$theta
    }

    list(
      sample = as.integer(sample_id),
      sourceSample = NA_integer_,
      source = "covariance",
      hypothesis = "simulation",
      modelLabel = "simulation",
      role = "simulation",
      theta = theta,
      omega = drawn$omega,
      sigma = base_params$sigma
    )
  })

  list(
    records = records,
    info = list(
      mode = "covariance",
      covarianceDraw = draw_mode,
      parameterPartition = list(
        drawn = unname(c(aligned$drawNames, drawn_omega_cols)),
        fixed = unname(c(
          setdiff(schema$thetaCols, aligned$drawNames),
          setdiff(schema$omegaCols, drawn_omega_cols),
          schema$sigmaCols
        )),
        # Persist WHY each block was held. cli_inform() above only reaches the
        # console; without this the reason is lost the moment the run ends, and
        # a saved run cannot explain why its OMEGA did not vary.
        heldOmegaBlocks = lapply(coverage$held, function(h) {
          list(
            etas = rownames(omega0)[h$index],
            reason = h$reason
          )
        })
      ),
      estimationInitialValues = "reference_fit"
    )
  )
}

.resolveSimulationParameters <- function(
  fit,
  samples,
  control,
  outputDir,
  schema
) {
  if (identical(control$parameterSource, "fixed")) {
    base_params <- .paramSetFromFit(fit)
    records <- lapply(seq_len(samples), function(sample_id) {
      list(
        sample = as.integer(sample_id),
        sourceSample = NA_integer_,
        source = "fixed",
        hypothesis = "simulation",
        modelLabel = "simulation",
        role = "simulation",
        theta = base_params$theta,
        omega = base_params$omega,
        sigma = base_params$sigma
      )
    })
    return(list(
      records = records,
      info = list(
        mode = "fixed",
        parameterPartition = list(
          drawn = character(0),
          fixed = unname(.schemaParamCols(schema))
        ),
        estimationInitialValues = "reference_fit"
      )
    ))
  }

  if (identical(control$parameterSource, "rawres")) {
    return(.resolveRawresParameterSets(
      fit = fit,
      samples = samples,
      control = control
    ))
  }

  .resolveCovarianceParameterSets(
    fit = fit,
    samples = samples,
    outputDir = outputDir,
    schema = schema,
    control = control
  )
}

.isObservationRow <- function(data) {
  mdv <- if ("MDV" %in% names(data)) {
    suppressWarnings(as.numeric(data$MDV))
  } else {
    rep(0, nrow(data))
  }
  evid <- if ("EVID" %in% names(data)) {
    suppressWarnings(as.numeric(data$EVID))
  } else {
    rep(0, nrow(data))
  }
  (is.na(evid) | evid == 0) & (is.na(mdv) | mdv == 0)
}

.applySimulatedValues <- function(refData, solved, appendColumns = NULL) {
  row_col <- .simulationRowIdCol
  if (!row_col %in% names(solved)) {
    .abortSSE("Simulation output did not preserve the internal row identifier.")
  }
  solved <- solved[order(solved[[row_col]]), , drop = FALSE]

  out <- refData
  obs <- .isObservationRow(out)
  solved_obs <- .isObservationRow(solved)

  if (sum(obs) != sum(solved_obs)) {
    .abortSSE(
      "Simulation output did not align with the original observation rows."
    )
  }
  if (!"sim" %in% names(solved)) {
    .abortSSE("Simulation output did not contain the {.field sim} column.")
  }

  out$DV[obs] <- solved$sim[solved_obs]

  for (column in appendColumns %||% character(0)) {
    if (!column %in% names(solved)) {
      .abortSSE(
        "Simulation output did not contain requested {.arg appendColumns} entry {.val {column}}."
      )
    }
    if (!column %in% names(out)) {
      out[[column]] <- NA
    }
    out[[column]][obs] <- solved[[column]][solved_obs]
  }

  out[[row_col]] <- NULL
  out
}

.callSimulationPostProcess <- function(
  postProcess,
  data,
  sample,
  paramSet,
  solved,
  referenceData,
  outputDir
) {
  if (is.null(postProcess)) {
    return(data)
  }

  available <- list(
    data = data,
    sample = as.integer(sample),
    paramSet = paramSet,
    solved = solved,
    referenceData = referenceData,
    outputDir = outputDir
  )
  fn_formals <- names(formals(postProcess))
  args <- if ("..." %in% fn_formals) {
    available
  } else {
    available[intersect(fn_formals, names(available))]
  }

  out <- tryCatch(
    do.call(postProcess, args),
    error = function(e) {
      .abortSSE(
        "simulationPostProcess failed for sample {.val {sample}}: {conditionMessage(e)}"
      )
    }
  )

  if (!is.data.frame(out)) {
    .abortSSE(
      "{.arg simulationPostProcess} must return a data frame for sample {.val {sample}}."
    )
  }

  as.data.frame(out, stringsAsFactors = FALSE)
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

.writeDatasetArtifacts <- function(data, dir, sample) {
  stem <- file.path(dir, sprintf("sim_%04d", sample))
  utils::write.csv(data, paste0(stem, ".csv"), row.names = FALSE, na = "NA")
  saveRDS(data, paste0(stem, ".rds"))
}

.fitArtifactPath <- function(dir, label, sample) {
  file.path(dir, sprintf("%s_%04d.rds", label, sample))
}

.fitErrorArtifact <- function(sample, label, error, call = NULL) {
  structure(
    list(
      sample = as.integer(sample),
      model_label = label,
      condition = structure(
        list(
          message = conditionMessage(error),
          class = class(error)
        ),
        class = "nlmixr2SSECondition"
      ),
      call = call,
      traceback = utils::capture.output(sys.calls())
    ),
    class = "nlmixr2SSEFitError"
  )
}

.omegaMatrixFromUi <- function(ui) {
  ini_df <- ui[["iniDf"]]
  if (!is.data.frame(ini_df)) {
    return(matrix(numeric(0), 0L, 0L))
  }

  omega_rows <- ini_df[
    !is.na(ini_df$neta1) & !is.na(ini_df$neta2),
    ,
    drop = FALSE
  ]
  if (nrow(omega_rows) == 0L) {
    return(matrix(numeric(0), 0L, 0L))
  }

  n_eta <- max(
    c(as.integer(omega_rows$neta1), as.integer(omega_rows$neta2)),
    na.rm = TRUE
  )
  mat <- matrix(0, n_eta, n_eta)
  diag_rows <- omega_rows[omega_rows$neta1 == omega_rows$neta2, , drop = FALSE]
  eta_names <- as.character(seq_len(n_eta))
  if (nrow(diag_rows) > 0L) {
    eta_names[diag_rows$neta1] <- diag_rows$name
  }
  dimnames(mat) <- list(eta_names, eta_names)

  est_vals <- if ("est" %in% names(omega_rows)) {
    as.numeric(omega_rows$est)
  } else {
    rep(0, nrow(omega_rows))
  }
  for (i in seq_len(nrow(omega_rows))) {
    row_idx <- as.integer(omega_rows$neta1[[i]])
    col_idx <- as.integer(omega_rows$neta2[[i]])
    mat[row_idx, col_idx] <- est_vals[[i]]
    mat[col_idx, row_idx] <- est_vals[[i]]
  }
  mat
}

.uiOmegaInfo <- function(ui) {
  ini_df <- ui[["iniDf"]]
  if (!is.data.frame(ini_df)) {
    return(data.frame(
      row = integer(0),
      col = integer(0),
      rowName = character(0),
      colName = character(0),
      fix = logical(0),
      stringsAsFactors = FALSE
    ))
  }

  omega_rows <- ini_df[
    !is.na(ini_df$neta1) & !is.na(ini_df$neta2) & ini_df$neta1 >= ini_df$neta2,
    ,
    drop = FALSE
  ]
  if (nrow(omega_rows) == 0L) {
    return(data.frame(
      row = integer(0),
      col = integer(0),
      rowName = character(0),
      colName = character(0),
      fix = logical(0),
      stringsAsFactors = FALSE
    ))
  }

  omega_rows <- omega_rows[
    order(omega_rows$neta2, omega_rows$neta1),
    ,
    drop = FALSE
  ]
  omega_mat <- .omegaMatrixFromUi(ui)

  data.frame(
    row = as.integer(omega_rows$neta1),
    col = as.integer(omega_rows$neta2),
    rowName = rownames(omega_mat)[as.integer(omega_rows$neta1)],
    colName = colnames(omega_mat)[as.integer(omega_rows$neta2)],
    fix = if ("fix" %in% names(omega_rows)) {
      !is.na(omega_rows$fix) & omega_rows$fix
    } else {
      FALSE
    },
    stringsAsFactors = FALSE
  )
}

.omegaValueForUi <- function(paramOmega, rowName, colName, row, col) {
  if (
    !is.matrix(paramOmega) || nrow(paramOmega) == 0L || ncol(paramOmega) == 0L
  ) {
    return(NA_real_)
  }

  if (!is.null(rownames(paramOmega)) && !is.null(colnames(paramOmega))) {
    if (
      rowName %in% rownames(paramOmega) && colName %in% colnames(paramOmega)
    ) {
      return(as.numeric(paramOmega[rowName, colName]))
    }
  }

  if (nrow(paramOmega) >= row && ncol(paramOmega) >= col) {
    return(as.numeric(paramOmega[row, col]))
  }

  NA_real_
}

.omegaRowExpressions <- function(current, updated) {
  if (!is.matrix(updated) || nrow(updated) == 0L || ncol(updated) == 0L) {
    return(list())
  }

  exprs <- list()
  for (i in seq_len(nrow(updated))) {
    current_vals <- as.numeric(current[i, seq_len(i), drop = TRUE])
    updated_vals <- as.numeric(updated[i, seq_len(i), drop = TRUE])
    if (identical(current_vals, updated_vals)) {
      next
    }

    rhs <- if (i == 1L) {
      updated_vals[[1L]]
    } else {
      as.call(c(list(as.name("c")), as.list(updated_vals)))
    }
    exprs[[length(exprs) + 1L]] <- call(
      "~",
      as.name(rownames(updated)[[i]]),
      rhs
    )
  }
  exprs
}

.uiWithParameterUpdates <- function(
  ui,
  paramSet,
  updateEstimated = FALSE,
  updateFixed = FALSE
) {
  if (
    !inherits(ui, "rxUi") ||
      (!isTRUE(updateEstimated) && !isTRUE(updateFixed))
  ) {
    return(ui)
  }

  ini_df <- ui[["iniDf"]]
  if (!is.data.frame(ini_df) || nrow(ini_df) == 0L) {
    return(ui)
  }

  exprs <- list()
  theta_rows <- which(
    !is.na(ini_df$ntheta) &
      !is.na(ini_df$name) &
      nzchar(ini_df$name)
  )
  theta_vals <- paramSet$theta %||% stats::setNames(numeric(0), character(0))

  for (idx in theta_rows) {
    name <- ini_df$name[[idx]]
    if (!name %in% names(theta_vals)) {
      next
    }
    fixed <- "fix" %in% names(ini_df) && isTRUE(ini_df$fix[[idx]])
    if ((fixed && isTRUE(updateFixed)) || (!fixed && isTRUE(updateEstimated))) {
      exprs[[length(exprs) + 1L]] <- call(
        "<-",
        as.name(name),
        theta_vals[[name]]
      )
    }
  }

  current_omega <- .omegaMatrixFromUi(ui)
  omega_info <- .uiOmegaInfo(ui)
  if (nrow(omega_info) > 0L) {
    updated_omega <- current_omega
    for (i in seq_len(nrow(omega_info))) {
      fixed <- isTRUE(omega_info$fix[[i]])
      if (!fixed && !isTRUE(updateEstimated)) {
        next
      }
      if (fixed && !isTRUE(updateFixed)) {
        next
      }

      value <- .omegaValueForUi(
        paramSet$omega,
        rowName = omega_info$rowName[[i]],
        colName = omega_info$colName[[i]],
        row = omega_info$row[[i]],
        col = omega_info$col[[i]]
      )
      if (is.na(value)) {
        next
      }
      updated_omega[
        omega_info$row[[i]],
        omega_info$col[[i]]
      ] <- value
      updated_omega[
        omega_info$col[[i]],
        omega_info$row[[i]]
      ] <- value
    }
    exprs <- c(exprs, .omegaRowExpressions(current_omega, updated_omega))
  }

  if (length(exprs) == 0L) {
    return(ui)
  }

  ui_obj <- ui
  eval_env <- new.env(parent = parent.frame())
  eval_env$ui_obj <- ui_obj
  eval_env$ini_eval_env <- eval_env
  if (requireNamespace("lotri", quietly = TRUE)) {
    eval_env$lotri <- lotri::lotri
  }
  call <- as.call(c(
    list(quote(rxode2::ini), quote(ui_obj)),
    exprs,
    list(envir = quote(ini_eval_env))
  ))
  suppressWarnings(suppressMessages(eval(call, envir = eval_env)))
}

.simulationRecord <- function(
  sample,
  fit,
  refData,
  control,
  outputDir,
  schema,
  paramSet
) {
  .ensureRxodeCache()

  sim_data <- refData
  sim_data[[.simulationRowIdCol]] <- seq_len(nrow(sim_data))

  solved <- nlmixr2utils::withRunSeed(
    outputDir,
    key = paste0("simulation-", sample),
    prefix = "sse",
    expr = {
      keep_cols <- unique(c(.simulationRowIdCol, control$appendColumns %||% character(0)))
      suppressWarnings(suppressMessages(local({
        utils::capture.output(
          result <- as.data.frame(
            rxode2::rxSolve(
              .getFitUi(fit),
              params = paramSet$theta,
              events = sim_data,
              omega = paramSet$omega,
              sigma = paramSet$sigma,
              addCov = TRUE,
              keep = keep_cols,
              returnType = "data.frame"
            ),
            stringsAsFactors = FALSE
          ),
          type = "output"
        )
        result
      })))
    }
  )

  sim_estimation_data <- .applySimulatedValues(
    refData = sim_data,
    solved = solved,
    appendColumns = control$appendColumns
  )
  sim_estimation_data <- .callSimulationPostProcess(
    postProcess = control$simulationPostProcess,
    data = sim_estimation_data,
    sample = sample,
    paramSet = paramSet,
    solved = solved,
    referenceData = refData,
    outputDir = outputDir
  )

  sim_param <- .simulationParamRecord(paramSet, schema)
  names(sim_param) <- .schemaParamCols(schema)

  if (isTRUE(control$saveDatasets)) {
    .writeDatasetArtifacts(
      sim_estimation_data,
      file.path(outputDir, "simulations"),
      sample
    )
  }

  list(
    sample = as.integer(sample),
    data = sim_estimation_data,
    initial = sim_param,
    paramSet = paramSet
  )
}

.fitTaskKey <- function(label, sample) {
  paste(label, sample, sep = "::")
}

.fitTaskRecord <- function(
  sample,
  spec,
  simRecord,
  unionSchema,
  outputDir,
  saveFits,
  control
) {
  .ensureRxodeCache()

  quiet_control <- nlmixr2utils::setQuietFastControl(spec$control)
  artifact_path <- .fitArtifactPath(
    file.path(outputDir, "fits"),
    spec$label,
    sample
  )
  fit_ui <- .uiWithParameterUpdates(
    spec$ui,
    simRecord$paramSet,
    updateEstimated = identical(control$parameterSource, "rawres") &&
      isTRUE(control$randomEstimationInits),
    updateFixed = identical(control$parameterSource, "rawres") &&
      isTRUE(control$updateFix)
  )

  fit_call <- quote(
    nlmixr2est::nlmixr2(
      fit_ui,
      simRecord$data,
      est = spec$est,
      control = quiet_control
    )
  )

  fit_result <- tryCatch(
    {
      fit_obj <- suppressWarnings(suppressMessages(local({
        utils::capture.output(result <- eval(fit_call), type = "output")
        result
      })))
      if (isTRUE(saveFits)) {
        saveRDS(fit_obj, artifact_path)
      }
      list(
        success = TRUE,
        row = nlmixr2utils::rawResultsRow(
          fit = fit_obj,
          source = "sse",
          hypothesis = spec$hypothesis,
          sample = sample,
          modelLabel = spec$label,
          role = spec$role,
          schema = unionSchema
        )
      )
    },
    error = function(e) {
      if (isTRUE(saveFits)) {
        saveRDS(
          .fitErrorArtifact(sample, spec$label, e, call = fit_call),
          artifact_path
        )
      }
      list(
        success = FALSE,
        row = nlmixr2utils::rawResultsRow(
          fit = NULL,
          source = "sse",
          hypothesis = spec$hypothesis,
          sample = sample,
          modelLabel = spec$label,
          role = spec$role,
          errorMessage = .truncateSSEMessage(conditionMessage(e)),
          objf = NA_real_,
          minimizationSuccessful = NA_integer_,
          covarianceStepSuccessful = NA_integer_,
          estimateNearBoundary = NA_integer_,
          significantDigits = NA_real_,
          conditionNumber = NA_real_,
          schema = unionSchema
        )
      )
    }
  )

  fit_result$sample <- as.integer(sample)
  fit_result$model_label <- spec$label
  fit_result$role <- spec$role
  fit_result$hypothesis <- spec$hypothesis
  fit_result
}

.initialEstimatesWide <- function(simRecords, schema, samples) {
  param_cols <- .schemaParamCols(schema)
  out <- data.frame(
    sample = as.integer(samples),
    matrix(
      NA_real_,
      nrow = length(samples),
      ncol = length(param_cols),
      dimnames = list(NULL, param_cols)
    ),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  for (record in simRecords) {
    idx <- match(record$sample, out$sample)
    out[idx, param_cols] <- as.numeric(record$initial[param_cols])
  }

  out
}

.wideToLongValues <- function(initialWide) {
  if (nrow(initialWide) == 0L || ncol(initialWide) <= 1L) {
    return(.emptyInitialValues())
  }
  param_cols <- setdiff(names(initialWide), "sample")
  values <- vector("list", length(param_cols))
  for (i in seq_along(param_cols)) {
    values[[i]] <- data.frame(
      replicate = as.integer(initialWide$sample),
      parameter = param_cols[[i]],
      value = as.numeric(initialWide[[param_cols[[i]]]]),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, values)
}

.wideToReferenceValues <- function(initialWide) {
  out <- .wideToLongValues(initialWide)
  names(out)[names(out) == "replicate"] <- "replicate"
  out[, c("parameter", "replicate", "value"), drop = FALSE]
}

.summaryMask <- function(rawResults, outFilter = NULL) {
  mask <- is.na(rawResults$error_message)
  if (!is.null(outFilter)) {
    predicate <- if (is.function(outFilter)) {
      outFilter
    } else {
      nlmixr2utils::setupRawResultsFilter(outFilter)
    }
    mask <- mask & predicate(rawResults)
  }
  mask[is.na(mask)] <- FALSE
  mask
}

.parameterStatsOrder <- c(
  "mean",
  "median",
  "sd",
  "max",
  "min",
  "skewness",
  "kurtosis",
  "rmse",
  "relative_rmse",
  "bias",
  "relative_bias",
  "relative_absolute_bias",
  "rse",
  "ci_0.5",
  "ci_2.5",
  "ci_5",
  "ci_95",
  "ci_97.5",
  "ci_99.5"
)

.ciProbMap <- c(
  "ci_0.5" = 0.005,
  "ci_2.5" = 0.025,
  "ci_5" = 0.05,
  "ci_95" = 0.95,
  "ci_97.5" = 0.975,
  "ci_99.5" = 0.995
)

.safeMean <- function(x) {
  if (length(x) == 0L) NA_real_ else mean(x)
}

.safeMedian <- function(x) {
  if (length(x) == 0L) NA_real_ else stats::median(x)
}

.safeSd <- function(x) {
  if (length(x) < 2L) NA_real_ else stats::sd(x)
}

.safeMin <- function(x) {
  if (length(x) == 0L) NA_real_ else min(x)
}

.safeMax <- function(x) {
  if (length(x) == 0L) NA_real_ else max(x)
}

.safeSkewness <- function(x) {
  n <- length(x)
  if (n < 3L) {
    return(NA_real_)
  }
  s <- stats::sd(x)
  if (!is.finite(s) || s == 0) {
    return(NA_real_)
  }
  m <- mean(x)
  (n / ((n - 1L) * (n - 2L))) * s^(-3) * sum((x - m)^3)
}

.safeKurtosis <- function(x) {
  n <- length(x)
  if (n < 4L) {
    return(NA_real_)
  }
  s <- stats::sd(x)
  if (!is.finite(s) || s == 0) {
    return(NA_real_)
  }
  m <- mean(x)
  (n * (n + 1L) / ((n - 1L) * (n - 2L) * (n - 3L))) *
    s^(-4) *
    sum((x - m)^4) -
    3 * (n - 1L)^2 / ((n - 2L) * (n - 3L))
}

.singleParameterSummary <- function(estimates, truth, matched) {
  stats <- stats::setNames(
    rep(NA_real_, length(.parameterStatsOrder)),
    .parameterStatsOrder
  )
  ns <- stats::setNames(
    rep(0L, length(.parameterStatsOrder)),
    .parameterStatsOrder
  )

  if (!isTRUE(matched)) {
    return(list(values = stats, n = ns))
  }

  estimate_ok <- !is.na(estimates)
  est <- estimates[estimate_ok]

  basic_n <- length(est)
  if (basic_n > 0L) {
    stats["mean"] <- .safeMean(est)
    stats["median"] <- .safeMedian(est)
    stats["sd"] <- .safeSd(est)
    stats["max"] <- .safeMax(est)
    stats["min"] <- .safeMin(est)
    stats["skewness"] <- .safeSkewness(est)
    stats["kurtosis"] <- .safeKurtosis(est)
    ns[c(
      "mean",
      "median",
      "sd",
      "max",
      "min",
      "skewness",
      "kurtosis"
    )] <- basic_n
  }

  pair_ok <- estimate_ok & !is.na(truth)
  est_pair <- estimates[pair_ok]
  true_pair <- truth[pair_ok]
  pair_n <- length(est_pair)

  if (pair_n > 0L) {
    diff <- est_pair - true_pair
    stats["rmse"] <- sqrt(mean(diff^2))
    stats["bias"] <- mean(diff)
    ns[c("rmse", "bias")] <- pair_n
  }

  rel_ok <- pair_ok & truth != 0
  est_rel <- estimates[rel_ok]
  true_rel <- truth[rel_ok]
  rel_n <- length(est_rel)

  if (rel_n > 0L) {
    rel_diff <- (est_rel - true_rel) / true_rel
    stats["relative_rmse"] <- 100 * sqrt(mean(rel_diff^2))
    stats["relative_bias"] <- 100 * mean(rel_diff)
    stats["relative_absolute_bias"] <- 100 * mean(abs(rel_diff))
    ns[c("relative_rmse", "relative_bias", "relative_absolute_bias")] <- rel_n
  }

  if (
    pair_n > 0L && all(!is.na(true_pair)) && length(unique(true_pair)) == 1L
  ) {
    true0 <- true_pair[[1L]]
    if (is.finite(true0) && true0 != 0) {
      stats["rse"] <- 100 * .safeSd(est_pair) / (true0 * sqrt(pair_n))
      ns["rse"] <- pair_n
      if (is.finite(stats["relative_bias"]) && is.finite(stats["rse"])) {
        for (nm in names(.ciProbMap)) {
          stats[[nm]] <- stats["relative_bias"] +
            stats["rse"] * stats::qnorm(.ciProbMap[[nm]])
          ns[[nm]] <- pair_n
        }
      }
    }
  }

  list(values = stats, n = ns)
}

.computeParameterSummary <- function(
  rawResults,
  initialWide,
  specMap,
  outFilter = NULL
) {
  if (nrow(rawResults) == 0L) {
    return(.emptyParameterSummary())
  }

  header <- attr(rawResults, "rawResultsHeader", exact = TRUE)
  param_cols <- c(header$theta_cols, header$omega_cols, header$sigma_cols)
  samples <- initialWide$sample
  truth_map <- initialWide
  row_mask <- .summaryMask(rawResults, outFilter)

  rows <- vector(
    "list",
    length(specMap) * length(param_cols) * length(.parameterStatsOrder)
  )
  row_index <- 1L

  for (label in names(specMap)) {
    spec_rows <- rawResults$sample > 0L &
      rawResults$model_label == label &
      row_mask
    role <- if (any(rawResults$model_label == label)) {
      rawResults$role[match(label, rawResults$model_label)]
    } else {
      NA_character_
    }
    model_params <- names(.schemaParamPresence(specMap[[label]]))

    for (param in param_cols) {
      matched <- param %in% model_params
      est <- if (matched && any(spec_rows)) {
        as.numeric(rawResults[spec_rows, param, drop = TRUE])
      } else {
        numeric(0)
      }
      sam <- if (matched && any(spec_rows)) {
        rawResults$sample[spec_rows]
      } else {
        integer(0)
      }
      truth <- if (matched && param %in% names(truth_map)) {
        truth_map[[param]][match(sam, samples)]
      } else {
        rep(NA_real_, length(est))
      }

      summary_one <- .singleParameterSummary(est, truth, matched = matched)

      for (stat_name in .parameterStatsOrder) {
        rows[[row_index]] <- data.frame(
          model_label = label,
          role = role,
          parameter = param,
          statistic = stat_name,
          value = as.numeric(summary_one$values[[stat_name]]),
          n_effective = as.integer(summary_one$n[[stat_name]]),
          matched = matched,
          stringsAsFactors = FALSE
        )
        row_index <- row_index + 1L
      }
    }
  }

  if (row_index == 1L) {
    return(.emptyParameterSummary())
  }
  do.call(rbind, rows[seq_len(row_index - 1L)])
}

.ofvThresholds <- c(0, 3.84, 5.99, 7.81, 9.49, 10.83)

.computeOfvSummary <- function(rawResults, fitSpecs, control, initialWide) {
  alternative_specs <- Filter(
    function(x) identical(x$role, "alternative"),
    fitSpecs
  )
  if (length(alternative_specs) == 0L || nrow(rawResults) == 0L) {
    return(.emptyOfvSummary())
  }

  samples <- initialWide$sample
  row_mask <- .summaryMask(rawResults, control$outFilter)
  filtered <- rawResults[row_mask, , drop = FALSE]

  ofv_by_label <- lapply(fitSpecs, function(spec) {
    rows <- filtered$sample > 0L & filtered$model_label == spec$label
    vec <- rep(NA_real_, length(samples))
    if (any(rows)) {
      vec[match(filtered$sample[rows], samples)] <- filtered$objf[rows]
    }
    vec
  })
  names(ofv_by_label) <- vapply(fitSpecs, `[[`, character(1), "label")

  ref_ofv <- if (!is.null(control$refOfv)) {
    rep(as.numeric(control$refOfv), length(samples))
  } else if (isTRUE(control$estimateSimulation)) {
    sim_label <- fitSpecs[[match(
      "simulation",
      vapply(
        fitSpecs,
        `[[`,
        character(1),
        "role"
      )
    )]]$label
    ofv_by_label[[sim_label]]
  } else {
    ofv_by_label[[alternative_specs[[1L]]$label]]
  }

  rows <- vector(
    "list",
    length(alternative_specs) * (1L + 2L * length(.ofvThresholds))
  )
  row_index <- 1L

  for (spec in alternative_specs) {
    delta <- ref_ofv - ofv_by_label[[spec$label]]
    valid <- is.finite(delta)

    rows[[row_index]] <- data.frame(
      model_label = spec$label,
      statistic = "mean_delta_ofv",
      threshold = NA_real_,
      direction = NA_character_,
      value = if (any(valid)) mean(delta[valid]) else NA_real_,
      stringsAsFactors = FALSE
    )
    row_index <- row_index + 1L

    for (threshold in .ofvThresholds) {
      rows[[row_index]] <- data.frame(
        model_label = spec$label,
        statistic = "pct_delta_above",
        threshold = threshold,
        direction = "type1",
        value = if (any(valid)) {
          100 * mean(delta[valid] > threshold)
        } else {
          NA_real_
        },
        stringsAsFactors = FALSE
      )
      row_index <- row_index + 1L

      rows[[row_index]] <- data.frame(
        model_label = spec$label,
        statistic = "pct_delta_below",
        threshold = threshold,
        direction = "power",
        value = if (any(valid)) {
          100 * mean(delta[valid] < (-threshold))
        } else {
          NA_real_
        },
        stringsAsFactors = FALSE
      )
      row_index <- row_index + 1L
    }
  }

  do.call(rbind, rows[seq_len(row_index - 1L)])
}

.wideParameterSummary <- function(parameterSummary) {
  if (is.null(parameterSummary) || nrow(parameterSummary) == 0L) {
    return(data.frame(
      model_label = character(0),
      role = character(0),
      parameter = character(0),
      matched = logical(0),
      stringsAsFactors = FALSE
    ))
  }

  wide <- stats::reshape(
    parameterSummary[, c(
      "model_label",
      "role",
      "parameter",
      "matched",
      "statistic",
      "value"
    )],
    idvar = c("model_label", "role", "parameter", "matched"),
    timevar = "statistic",
    direction = "wide"
  )

  names(wide) <- sub("^value\\.", "", names(wide))
  ordered_stats <- intersect(.parameterStatsOrder, names(wide))
  wide <- wide[,
    c("model_label", "role", "parameter", "matched", ordered_stats),
    drop = FALSE
  ]
  wide[order(wide$model_label, wide$role, wide$parameter), , drop = FALSE]
}

.wideOfvSummary <- function(ofvSummary) {
  if (nrow(ofvSummary) == 0L) {
    return(data.frame(model_label = character(0), stringsAsFactors = FALSE))
  }
  ofvSummary <- ofvSummary[!is.na(ofvSummary$statistic), , drop = FALSE]
  if (nrow(ofvSummary) == 0L) {
    return(data.frame(model_label = character(0), stringsAsFactors = FALSE))
  }

  metric <- ifelse(
    is.na(ofvSummary$threshold),
    ofvSummary$statistic,
    paste0(ofvSummary$statistic, "_", ofvSummary$threshold)
  )
  wide <- stats::reshape(
    data.frame(
      model_label = ofvSummary$model_label,
      metric = metric,
      value = ofvSummary$value,
      stringsAsFactors = FALSE
    ),
    idvar = "model_label",
    timevar = "metric",
    direction = "wide"
  )
  names(wide) <- sub("^value\\.", "", names(wide))

  desired <- c(
    "mean_delta_ofv",
    paste0("pct_delta_above_", .ofvThresholds),
    paste0("pct_delta_below_", .ofvThresholds)
  )
  desired <- intersect(desired, names(wide))
  wide <- wide[, c("model_label", desired), drop = FALSE]
  wide[order(wide$model_label), , drop = FALSE]
}

.writeInitialEstimates <- function(initialWide, dir) {
  utils::write.csv(
    initialWide,
    file.path(dir, "initial_estimates.csv"),
    row.names = FALSE,
    na = "NA"
  )
}

.writeSseResults <- function(parameterSummary, dir, basename = "sse_results") {
  utils::write.csv(
    .wideParameterSummary(parameterSummary),
    file.path(dir, paste0(basename, ".csv")),
    row.names = FALSE,
    na = "NA"
  )
}

.writeSseSummary <- function(
  runInfo,
  referenceValues,
  initialValues,
  parameterSummary,
  ofvSummary,
  dir,
  basename = "sse_summary"
) {
  saveRDS(
    list(
      runInfo = runInfo,
      referenceValues = referenceValues,
      initialValues = initialValues,
      parameterSummary = parameterSummary,
      ofvSummary = ofvSummary,
      powerSummary = ofvSummary[!is.na(ofvSummary$direction) & ofvSummary$direction == "power", , drop = FALSE]
    ),
    file.path(dir, paste0(basename, ".rds"))
  )
}

.readInitialEstimates <- function(dir) {
  path <- file.path(dir, "initial_estimates.csv")
  if (!file.exists(path)) {
    return(data.frame(sample = integer(0), stringsAsFactors = FALSE))
  }
  utils::read.csv(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = "NA"
  )
}

.loadCompletedSSE <- function(outputDir, runInfo, alternativeSpecs = NULL) {
  raw_results <- nlmixr2utils::readRawResults(outputDir)
  initial_wide <- .readInitialEstimates(outputDir)
  if (is.null(alternativeSpecs)) {
    alternativeSpecs <- .alternativeSnapshotsFromRunInfo(runInfo)
  }
  summary_path <- file.path(outputDir, "sse_summary.rds")
  summary_payload <- if (file.exists(summary_path)) {
    readRDS(summary_path)
  } else {
    list(
      referenceValues = .wideToReferenceValues(initial_wide),
      initialValues = .wideToLongValues(initial_wide),
      parameterSummary = .emptyParameterSummary(),
      ofvSummary = .emptyOfvSummary(),
      powerSummary = .emptyOfvSummary()
    )
  }

  .newNlmixr2SSE(
    runInfo = runInfo,
    rawResults = raw_results,
    alternativeSpecs = alternativeSpecs,
    outputDir = outputDir,
    timestamp = Sys.time(),
    referenceValues = summary_payload$referenceValues %||% .emptyReferenceValues(),
    initialValues = summary_payload$initialValues %||% .emptyInitialValues(),
    parameterSummary = summary_payload$parameterSummary %||% .emptyParameterSummary(),
    ofvSummary = summary_payload$ofvSummary %||% .emptyOfvSummary(),
    powerSummary = summary_payload$powerSummary %||% .emptyOfvSummary()
  )
}

.newNlmixr2SSE <- function(
  runInfo,
  rawResults,
  alternativeSpecs,
  outputDir,
  timestamp,
  referenceValues = .emptyReferenceValues(),
  initialValues = .emptyInitialValues(),
  parameterSummary = .emptyParameterSummary(),
  ofvSummary = .emptyOfvSummary(),
  powerSummary = .emptyOfvSummary()
) {
  structure(
    list(
      runInfo = runInfo,
      referenceValues = referenceValues,
      initialValues = initialValues,
      rawResults = rawResults,
      parameterSummary = parameterSummary,
      ofvSummary = ofvSummary,
      powerSummary = powerSummary,
      alternativeSpecs = alternativeSpecs,
      outputDir = outputDir,
      timestamp = timestamp
    ),
    class = c("nlmixr2SSE", "list")
  )
}
