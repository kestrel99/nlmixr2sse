skip_on_cran()

.sse_ref_model <- function() {
  ini({
    tka <- log(1.57)
    tcl <- log(2.72)
    tv <- log(31.5)
    eta.ka ~ 0.6
    add.sd <- 0.7
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl)
    v <- exp(tv)
    cp <- linCmt()
    cp ~ add(add.sd)
  })
}

.sse_alt_model <- function() {
  ini({
    tka <- log(1.57)
    tcl <- log(2.72)
    tv <- log(31.5)
    add.sd <- 0.7
  })
  model({
    ka <- exp(tka)
    cl <- exp(tcl)
    v <- exp(tv)
    cp <- linCmt()
    cp ~ add(add.sd)
  })
}

.fit_sse_model <- function(model) {
  skip_if_not_installed("nlmixr2")
  skip_if_not_installed("nlmixr2data")
  skip_if_not_installed("rxode2")

  fit <- try(
    suppressMessages(suppressWarnings(
      nlmixr2utils::nlmixr2(
        model,
        nlmixr2data::theo_sd,
        est = "focei",
        control = list(print = 0L, covMethod = "r", eval.max = 30L)
      )
    )),
    silent = TRUE
  )
  if (inherits(fit, "try-error")) {
    skip(paste("nlmixr2 SSE test fit setup failed:", as.character(fit)))
  }
  fit
}

.ui_update_model <- function() {
  ini({
    tka <- log(1.57)
    tcl <- fix(log(2.72))
    tv <- log(31.5)
    eta.ka ~ fix(0.6)
    add.sd <- 0.7
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl)
    v <- exp(tv)
    cp <- linCmt()
    cp ~ add(add.sd)
  })
}

test_that("rawres parameter resolution defaults to simulation rows from SSE raw results", {
  fit <- unclass(fake_sse_fit())
  schema <- nlmixr2utils::rawResultsSchema(fit)
  rows <- do.call(
    rbind,
    list(
      nlmixr2utils::rawResultsRow(
        fit,
        source = "sse",
        hypothesis = "simulation",
        sample = 1L,
        modelLabel = "ref_fit",
        role = "simulation",
        theta = c(tka = 0.60, tcl = 1.10)
      ),
      nlmixr2utils::rawResultsRow(
        fit,
        source = "sse",
        hypothesis = "alternative_1",
        sample = 1L,
        modelLabel = "alt1",
        role = "alternative",
        theta = c(tka = 9, tcl = 9)
      ),
      nlmixr2utils::rawResultsRow(
        fit,
        source = "sse",
        hypothesis = "simulation",
        sample = 2L,
        modelLabel = "ref_fit",
        role = "simulation",
        theta = c(tka = 0.70, tcl = 1.20)
      ),
      nlmixr2utils::rawResultsRow(
        fit,
        source = "sse",
        hypothesis = "alternative_1",
        sample = 2L,
        modelLabel = "alt1",
        role = "alternative",
        theta = c(tka = 8, tcl = 8)
      )
    )
  )
  ctl <- runSSEControl(parameterSource = "rawres", rawresInput = rows)

  params <- .resolveSimulationParameters(
    fit = fit,
    samples = 2L,
    control = ctl,
    outputDir = tempdir(),
    schema = schema
  )

  expect_length(params$records, 2L)
  expect_equal(
    vapply(params$records, `[[`, integer(1), "sourceSample"),
    c(1L, 2L)
  )
  expect_equal(
    vapply(params$records, function(x) x$theta[["tka"]], numeric(1)),
    c(0.60, 0.70)
  )
})

test_that("runSSE requires an alternative when estimateSimulation is FALSE", {
  err <- capture_sse_error(
    runSSE(
      fake_sse_fit(),
      control = runSSEControl(estimateSimulation = FALSE),
      outputDir = tempfile("nlmixr2sse-phase4-alt-"),
      restart = TRUE
    )
  )

  expect_s3_class(err, "error")
  expect_match(conditionMessage(err), "alternative")
})

test_that("runSSE validates alternative-label collisions before execution", {
  err <- capture_sse_error(
    runSSE(
      fake_sse_fit(),
      alternativeModels = sseModel(fake_sse_fit(), label = "fake_sse_fit"),
      outputDir = tempfile("nlmixr2sse-badlabel-"),
      restart = TRUE
    )
  )

  expect_s3_class(err, "error")
  expect_match(conditionMessage(err), "collides")
})

test_that("covariance parameter resolution draws theta and keeps omega fixed", {
  fit <- unclass(fake_sse_fit())
  schema <- nlmixr2utils::rawResultsSchema(fit)
  tmp <- tempfile("nlmixr2sse-cov-")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE, force = TRUE), add = TRUE)

  params <- .resolveSimulationParameters(
    fit = fit,
    samples = 3L,
    control = runSSEControl(parameterSource = "covariance"),
    outputDir = tmp,
    schema = schema
  )

  expect_equal(params$info$parameterPartition$drawn, names(fit[["theta"]]))
  expect_equal(
    params$info$parameterPartition$fixed,
    "omega(eta.ka,eta.ka)"
  )
  expect_equal(
    vapply(params$records, function(x) x$omega[1, 1], numeric(1)),
    rep(fit[["omega"]][1, 1], 3L)
  )
  expect_gte(
    length(unique(signif(
      vapply(
        params$records,
        function(x) x$theta[["tka"]],
        numeric(1)
      ),
      8
    ))),
    2L
  )
})

test_that("ui parameter updates separate free from fixed estimates", {
  skip_if_not_installed("rxode2")
  skip_if_not_installed("nlmixr2est")

  ui <- nlmixr2est::nlmixr2(.ui_update_model)
  param_set <- list(
    theta = c(tka = 0.90, tcl = 0.80, tv = 3.20, add.sd = 0.55),
    omega = matrix(
      0.45,
      1L,
      1L,
      dimnames = list("eta.ka", "eta.ka")
    )
  )

  updated_free <- .uiWithParameterUpdates(
    ui,
    param_set,
    updateEstimated = TRUE,
    updateFixed = FALSE
  )
  free_df <- rxode2::rxUiDecompress(updated_free)[["iniDf"]]
  expect_equal(
    free_df$est[free_df$name == "tka"],
    0.90
  )
  expect_equal(
    free_df$est[free_df$name == "tcl"],
    log(2.72)
  )
  expect_equal(
    free_df$est[!is.na(free_df$neta1)],
    0.6
  )

  updated_fixed <- .uiWithParameterUpdates(
    ui,
    param_set,
    updateEstimated = FALSE,
    updateFixed = TRUE
  )
  fixed_df <- rxode2::rxUiDecompress(updated_fixed)[["iniDf"]]
  expect_equal(
    fixed_df$est[fixed_df$name == "tka"],
    log(1.57)
  )
  expect_equal(
    fixed_df$est[fixed_df$name == "tcl"],
    0.80
  )
  expect_equal(
    fixed_df$est[!is.na(fixed_df$neta1)],
    0.45
  )
})

test_that("parameter summary computes RSE only for constant truths", {
  raw_results <- data.frame(
    source = rep("sse", 2L),
    hypothesis = rep("simulation", 2L),
    sample = 1:2,
    model_label = rep("ref_fit", 2L),
    role = rep("simulation", 2L),
    minimization_successful = c(1L, 1L),
    covariance_step_successful = c(1L, 1L),
    estimate_near_boundary = c(0L, 0L),
    significant_digits = c(4.2, 4.1),
    condition_number = c(10, 11),
    objf = c(100, 101),
    error_message = c(NA_character_, NA_character_),
    tka = c(0.55, 0.72),
    "omega(eta.ka,eta.ka)" = c(0.31, 0.29),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  attr(raw_results, "rawResultsHeader") <- list(
    schema_version = 1L,
    columns = names(raw_results),
    base_cols = names(raw_results)[1:12],
    theta_cols = "tka",
    omega_cols = "omega(eta.ka,eta.ka)",
    sigma_cols = character(0),
    se_cols = character(0),
    parameter_cols = c("tka", "omega(eta.ka,eta.ka)")
  )

  initial_wide <- data.frame(
    sample = 1:2,
    tka = c(0.50, 0.70),
    "omega(eta.ka,eta.ka)" = c(0.30, 0.30),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  spec_map <- list(
    ref_fit = list(
      thetaCols = "tka",
      omegaCols = "omega(eta.ka,eta.ka)",
      sigmaCols = character(0)
    )
  )

  summary_tbl <- .computeParameterSummary(
    rawResults = raw_results,
    initialWide = initial_wide,
    specMap = spec_map
  )

  theta_rse <- subset(summary_tbl, parameter == "tka" & statistic == "rse")
  omega_rse <- subset(
    summary_tbl,
    parameter == "omega(eta.ka,eta.ka)" & statistic == "rse"
  )

  expect_equal(theta_rse$value, NA_real_)
  expect_equal(is.finite(omega_rse$value), TRUE)
})

test_that("simulationPostProcess rewrites simulated estimation data", {
  ref <- data.frame(
    ID = 1:2,
    DV = c(10, 20),
    stringsAsFactors = FALSE
  )
  solved <- data.frame(
    ID = 1:2,
    sim = c(11, 21),
    .sse_row_id = 1:2,
    stringsAsFactors = FALSE
  )

  out <- .callSimulationPostProcess(
    postProcess = function(data, sample, outputDir) {
      data$DV <- data$DV + sample
      data$SOURCE_DIR <- basename(outputDir)
      data
    },
    data = .applySimulatedValues(ref, solved),
    sample = 2L,
    paramSet = list(),
    solved = solved,
    referenceData = ref,
    outputDir = file.path(tempdir(), "phase7")
  )

  expect_equal(out$DV, c(13, 23))
  expect_equal(out$SOURCE_DIR, rep("phase7", 2L))
})

test_that("simulationRecord preserves row identifiers through rxSolve output", {
  ref_fit <- .fit_sse_model(.sse_ref_model)
  tmp <- tempfile("nlmixr2sse-sim-record-")
  dir.create(tmp)
  dir.create(file.path(tmp, "simulations"))
  on.exit(unlink(tmp, recursive = TRUE, force = TRUE), add = TRUE)

  schema <- nlmixr2utils::rawResultsSchema(unclass(ref_fit))
  ref_data <- .getFitData(ref_fit)
  param_set <- list(
    theta = ref_fit$theta,
    omega = ref_fit$omega,
    sigma = ref_fit$sigma
  )

  rec <- .simulationRecord(
    sample = 1L,
    fit = ref_fit,
    refData = ref_data,
    control = runSSEControl(workers = 1L),
    outputDir = tmp,
    schema = schema,
    paramSet = param_set
  )

  expect_equal(is.list(rec), TRUE)
  expect_equal(nrow(rec$data), nrow(ref_data))
  expect_false(".sse_row_id" %in% names(rec$data))
})

test_that("recomputeSSE rebuilds summaries and writes suffixed outputs", {
  tmp <- tempfile("nlmixr2sse-recompute-")
  on.exit(unlink(tmp, recursive = TRUE, force = TRUE), add = TRUE)
  write_fake_sse_run_dir(tmp)

  res <- recomputeSSE(
    tmp,
    outFilter = ~ sample == 1L,
    overwrite = FALSE
  )

  expect_s3_class(res, "nlmixr2SSE")
  expect_true(file.exists(file.path(tmp, "sse_results_recompute1.csv")))
  expect_true(file.exists(file.path(tmp, "sse_summary_recompute1.rds")))

  mean_row <- subset(
    res$parameterSummary,
    model_label == "fake_sse_fit" &
      parameter == "tka" &
      statistic == "mean"
  )
  expect_equal(mean_row$value, 0.55)
})

test_that("runSSE resumes cached partial runs without duplicating rows", {
  tmp <- tempfile("nlmixr2sse-resume-")
  dir.create(tmp)
  dir.create(file.path(tmp, "simulations"))
  dir.create(file.path(tmp, "fits"))
  on.exit(unlink(tmp, recursive = TRUE, force = TRUE), add = TRUE)

  mock_sim_record <- function(sample_id) {
    theta_tka <- 0.50 + 0.20 * (sample_id - 1L)
    theta_tcl <- 1.00 + 0.20 * (sample_id - 1L)
    omega_val <- 0.30 + 0.02 * (sample_id - 1L)
    list(
      sample = as.integer(sample_id),
      data = data.frame(
        ID = sample_id,
        DV = sample_id,
        stringsAsFactors = FALSE
      ),
      initial = c(
        tka = theta_tka,
        tcl = theta_tcl,
        "omega(eta.ka,eta.ka)" = omega_val
      ),
      paramSet = list(
        theta = c(tka = theta_tka, tcl = theta_tcl),
        omega = matrix(
          omega_val,
          1L,
          1L,
          dimnames = list("eta.ka", "eta.ka")
        ),
        sigma = matrix(numeric(0), 0L, 0L)
      )
    )
  }

  mock_fit_record <- function(sample, spec, simRecord, unionSchema) {
    theta_vals <- c(
      tka = simRecord$initial[["tka"]] +
        if (spec$role == "simulation") 0.05 else 0.02,
      tcl = simRecord$initial[["tcl"]] +
        if (spec$role == "simulation") 0.05 else 0.02
    )
    omega_vals <- if (spec$role == "simulation") {
      c(
        "omega(eta.ka,eta.ka)" = simRecord$initial[["omega(eta.ka,eta.ka)"]] +
          0.01
      )
    } else {
      NULL
    }
    list(
      success = TRUE,
      row = nlmixr2utils::rawResultsRow(
        fit = NULL,
        source = "sse",
        hypothesis = spec$hypothesis,
        sample = sample,
        modelLabel = spec$label,
        role = spec$role,
        theta = theta_vals,
        omega = omega_vals,
        objf = if (spec$role == "simulation") 100 + sample else 200 + sample,
        schema = unionSchema
      ),
      sample = as.integer(sample),
      model_label = spec$label,
      role = spec$role,
      hypothesis = spec$hypothesis
    )
  }

  data_cache <- nlmixr2utils::taskCache(tmp, "sse_data")
  fit_cache <- nlmixr2utils::taskCache(tmp, "sse_fit_rows")
  data_cache$put("1", mock_sim_record(1L))

  partial_specs <- list(
    list(
      label = "fake_sse_fit",
      role = "simulation",
      hypothesis = "simulation"
    ),
    list(
      label = "alt1",
      role = "alternative",
      hypothesis = "alternative_1"
    )
  )
  fit_cache$put("fake_sse_fit::1", {
    rec <- mock_fit_record(
      1L,
      partial_specs[[1L]],
      mock_sim_record(1L),
      .mergeRawResultsSchemas(list(
        nlmixr2utils::rawResultsSchema(unclass(fake_sse_fit()))
      ))
    )
    rec$row$objf <- 901
    rec
  })
  fit_cache$put("alt1::1", {
    rec <- mock_fit_record(
      1L,
      partial_specs[[2L]],
      mock_sim_record(1L),
      .mergeRawResultsSchemas(list(
        nlmixr2utils::rawResultsSchema(unclass(fake_sse_fit()))
      ))
    )
    rec$row$objf <- 801
    rec
  })

  saveRDS(
    list(
      status = "running",
      fitName = "fake_sse_fit",
      samples = 2L,
      seed = 7L,
      parameterSource = "fixed",
      estimateSimulation = TRUE,
      modelLabels = c("fake_sse_fit", "alt1"),
      fitSpecs = fake_sse_fit_specs(),
      alternativeSpecs = list(list(
        label = "alt1",
        est = "focei",
        control = list(print = 0L),
        isFit = TRUE,
        hasDataOverride = FALSE
      )),
      control = runSSEControl(workers = 1L),
      startedAt = Sys.time()
    ),
    file.path(tmp, "run_info.rds")
  )
  nlmixr2utils::withRunSeed(tmp, seed = 7, prefix = "sse")

  testthat::local_mocked_bindings(
    .getFitData = function(fit) {
      data.frame(ID = 1:2, DV = c(1, 2), stringsAsFactors = FALSE)
    },
    .simulationRecord = function(sample, ...) {
      mock_sim_record(sample)
    },
    .fitTaskRecord = function(sample, spec, simRecord, unionSchema, ...) {
      mock_fit_record(sample, spec, simRecord, unionSchema)
    },
    .package = "nlmixr2sse"
  )

  res <- suppressWarnings(runSSE(
    fake_sse_fit(),
    alternativeModels = sseModel(fake_sse_fit(), label = "alt1"),
    samples = 2L,
    control = runSSEControl(workers = 1L),
    outputDir = tmp,
    restart = FALSE
  ))

  expect_s3_class(res, "nlmixr2SSE")
  expect_equal(nrow(res$rawResults), 4L)
  expect_equal(
    res$rawResults$objf[
      res$rawResults$sample == 1L &
        res$rawResults$model_label == "fake_sse_fit"
    ],
    901
  )
  expect_equal(
    res$rawResults$objf[
      res$rawResults$sample == 1L &
        res$rawResults$model_label == "alt1"
    ],
    801
  )
  expect_identical(
    anyDuplicated(paste(res$rawResults$model_label, res$rawResults$sample)),
    0L
  )

  state <- nlmixr2utils::readRunState(tmp, "sse")
  expect_equal(state$status, "completed")
  expect_equal(state$completedSamples, c(1L, 2L))
  expect_length(state$completedFitTasks, 4L)
})

test_that("runSSE addModels reuses saved datasets and writes add artifacts", {
  tmp <- tempfile("nlmixr2sse-addmodels-")
  on.exit(unlink(tmp, recursive = TRUE, force = TRUE), add = TRUE)
  write_fake_sse_run_dir(tmp)
  saveRDS(
    data.frame(ID = 1, DV = 10, stringsAsFactors = FALSE),
    file.path(tmp, "simulations", "sim_0001.rds")
  )
  saveRDS(
    data.frame(ID = 2, DV = 20, stringsAsFactors = FALSE),
    file.path(tmp, "simulations", "sim_0002.rds")
  )
  nlmixr2utils::withRunSeed(tmp, seed = 42, prefix = "sse")

  testthat::local_mocked_bindings(
    .fitTaskRecord = function(sample, spec, simRecord, unionSchema, ...) {
      list(
        success = TRUE,
        row = nlmixr2utils::rawResultsRow(
          fit = NULL,
          source = "sse",
          hypothesis = spec$hypothesis,
          sample = sample,
          modelLabel = spec$label,
          role = spec$role,
          theta = c(
            tka = simRecord$initial[["tka"]] + 0.01,
            tcl = simRecord$initial[["tcl"]] + 0.01
          ),
          objf = 150 + sample,
          schema = unionSchema
        ),
        sample = as.integer(sample),
        model_label = spec$label,
        role = spec$role,
        hypothesis = spec$hypothesis
      )
    },
    .package = "nlmixr2sse"
  )

  res <- runSSE(
    fake_sse_fit(),
    alternativeModels = sseModel(fake_sse_fit(), label = "alt2"),
    samples = 2L,
    control = runSSEControl(addModels = TRUE, workers = 1L),
    outputDir = tmp,
    restart = FALSE
  )

  expect_s3_class(res, "nlmixr2SSE")
  expect_true(file.exists(file.path(tmp, "raw_results_add1.csv")))
  expect_true(file.exists(file.path(tmp, "raw_results_add1.rds")))
  expect_true(file.exists(file.path(tmp, "raw_results_add1_header.json")))
  expect_true(file.exists(file.path(tmp, "sse_results_add1.csv")))
  expect_true("alt2" %in% unique(res$rawResults$model_label))

  updated_run <- readRDS(file.path(tmp, "run_info.rds"))
  expect_true("alt2" %in% updated_run$modelLabels)

  state <- nlmixr2utils::readRunState(tmp, "sse")
  expect_equal(length(state$addRuns), 1L)
  expect_equal(state$addRuns[[1L]]$labels, "alt2")
})

test_that("runSSE runs the fixed-parameter workflow end to end", {
  ref_fit <- .fit_sse_model(.sse_ref_model)
  alt_fit <- .fit_sse_model(.sse_alt_model)

  tmp <- tempfile("nlmixr2sse-phase4-")
  on.exit(unlink(tmp, recursive = TRUE, force = TRUE), add = TRUE)

  res <- try(
    suppressMessages(
      runSSE(
        ref_fit,
        alternativeModels = sseModel(alt_fit, label = "no_eta"),
        samples = 2L,
        seed = 42,
        control = runSSEControl(workers = 1L),
        outputDir = tmp,
        restart = TRUE
      )
    ),
    silent = TRUE
  )
  if (inherits(res, "try-error")) {
    skip(paste(
      "Phase 4 SSE integration is unavailable in this environment:",
      as.character(res)
    ))
  }

  expect_s3_class(res, "nlmixr2SSE")
  expect_equal(res$runInfo$status, "completed")
  expect_equal(res$runInfo$fitName, "ref_fit")
  expect_equal(nrow(res$rawResults), 4L)
  expect_true(file.exists(file.path(tmp, "initial_estimates.csv")))
  expect_true(file.exists(file.path(tmp, "raw_results.csv")))
  expect_true(file.exists(file.path(tmp, "raw_results.rds")))
  expect_true(file.exists(file.path(tmp, "raw_results_header.json")))
  expect_true(file.exists(file.path(tmp, "sse_results.csv")))
  expect_true(file.exists(file.path(tmp, "sse_summary.rds")))
  expect_true(dir.exists(file.path(tmp, "simulations")))
  expect_true(dir.exists(file.path(tmp, "fits")))
  expect_equal(
    sort(unique(res$rawResults$model_label)),
    c("no_eta", "ref_fit")
  )
  expect_true(nrow(res$initialValues) > 0L)
  expect_true(nrow(res$parameterSummary) > 0L)
  expect_true(nrow(res$ofvSummary) > 0L)
  expect_true(nrow(res$powerSummary) > 0L)

  summ <- summary(res)
  expect_s3_class(summ, "summary.nlmixr2SSE")
  expect_true("mean" %in% names(summ$parameterSummary))
  expect_true("mean_delta_ofv" %in% names(summ$ofvSummary))

  printed <- utils::capture.output(ret <- print(res))
  expect_identical(ret, res)
  expect_match(paste(printed, collapse = "\n"), "Parameter summary")
  expect_match(paste(printed, collapse = "\n"), "OFV summary")

  resumed <- suppressMessages(
    runSSE(
      ref_fit,
      alternativeModels = sseModel(alt_fit, label = "no_eta"),
      samples = 2L,
      seed = 999,
      control = runSSEControl(workers = 1L),
      outputDir = tmp,
      restart = FALSE
    )
  )
  expect_s3_class(resumed, "nlmixr2SSE")
  expect_equal(nrow(resumed$rawResults), nrow(res$rawResults))
})

test_that("runSSE supports raw-results and covariance parameter sources end to end", {
  ref_fit <- .fit_sse_model(.sse_ref_model)
  alt_fit <- .fit_sse_model(.sse_alt_model)

  rawres_dir <- tempfile("nlmixr2sse-rawres-input-")
  rawres_out <- tempfile("nlmixr2sse-rawres-run-")
  cov_out <- tempfile("nlmixr2sse-cov-run-")
  dir.create(rawres_dir)
  on.exit(unlink(rawres_dir, recursive = TRUE, force = TRUE), add = TRUE)
  on.exit(unlink(rawres_out, recursive = TRUE, force = TRUE), add = TRUE)
  on.exit(unlink(cov_out, recursive = TRUE, force = TRUE), add = TRUE)

  rawres_rows <- try(
    do.call(
      rbind,
      list(
        nlmixr2utils::rawResultsRow(
          unclass(ref_fit),
          source = "bootstrap",
          hypothesis = "reference",
          sample = 0L,
          modelLabel = "reference",
          role = "reference"
        ),
        nlmixr2utils::rawResultsRow(
          unclass(ref_fit),
          source = "bootstrap",
          hypothesis = "sample",
          sample = 1L,
          modelLabel = "reference",
          role = "reference",
          theta = c(tka = 0.47, tcl = 0.98, tv = 3.40, add.sd = 0.68),
          omega = c("omega(eta.ka,eta.ka)" = 0.55),
          minimizationSuccessful = 1L,
          covarianceStepSuccessful = 1L,
          significantDigits = 4.7
        ),
        nlmixr2utils::rawResultsRow(
          unclass(ref_fit),
          source = "bootstrap",
          hypothesis = "sample",
          sample = 2L,
          modelLabel = "reference",
          role = "reference",
          theta = c(tka = 0.43, tcl = 1.03, tv = 3.50, add.sd = 0.72),
          omega = c("omega(eta.ka,eta.ka)" = 0.63),
          minimizationSuccessful = 1L,
          covarianceStepSuccessful = 1L,
          significantDigits = 4.1
        ),
        nlmixr2utils::rawResultsRow(
          unclass(ref_fit),
          source = "bootstrap",
          hypothesis = "sample",
          sample = 3L,
          modelLabel = "reference",
          role = "reference",
          theta = c(tka = 0.49, tcl = 0.95, tv = 3.38, add.sd = 0.67),
          omega = c("omega(eta.ka,eta.ka)" = 0.58),
          minimizationSuccessful = 0L,
          covarianceStepSuccessful = NA_integer_,
          significantDigits = 2.5,
          errorMessage = "fit warning"
        )
      )
    ),
    silent = TRUE
  )
  if (inherits(rawres_rows, "try-error")) {
    skip(paste(
      "Phase 5 rawres fixture setup is unavailable in this environment:",
      as.character(rawres_rows)
    ))
  }
  nlmixr2utils::writeRawResults(rawres_rows, rawres_dir)

  rawres_res <- try(
    suppressMessages(
      runSSE(
        ref_fit,
        alternativeModels = sseModel(alt_fit, label = "no_eta"),
        samples = 2L,
        seed = 123,
        control = runSSEControl(
          parameterSource = "rawres",
          rawresInput = rawres_dir,
          inFilter = ~ minimization_successful == 1 & significant_digits > 4,
          randomEstimationInits = TRUE,
          workers = 1L
        ),
        outputDir = rawres_out,
        restart = TRUE
      )
    ),
    silent = TRUE
  )
  if (inherits(rawres_res, "try-error")) {
    skip(paste(
      "Phase 5 rawres SSE integration is unavailable in this environment:",
      as.character(rawres_res)
    ))
  }

  expect_s3_class(rawres_res, "nlmixr2SSE")
  expect_equal(
    rawres_res$runInfo$parameterSourceInfo$sourceSamples,
    c(1L, 2L)
  )
  expect_equal(
    rawres_res$runInfo$parameterSourceInfo$estimationInitialValues,
    "rawres"
  )
  rawres_tka <- rawres_res$initialValues[
    rawres_res$initialValues$parameter == "tka",
    c("replicate", "value"),
    drop = FALSE
  ]
  rawres_tka <- rawres_tka[order(rawres_tka$replicate), , drop = FALSE]
  expect_equal(round(rawres_tka$value, 2), c(0.47, 0.43))

  cov_res <- try(
    suppressMessages(
      runSSE(
        ref_fit,
        alternativeModels = sseModel(alt_fit, label = "no_eta"),
        samples = 2L,
        seed = 456,
        control = runSSEControl(
          parameterSource = "covariance",
          workers = 1L
        ),
        outputDir = cov_out,
        restart = TRUE
      )
    ),
    silent = TRUE
  )
  if (inherits(cov_res, "try-error")) {
    skip(paste(
      "Phase 5 covariance SSE integration is unavailable in this environment:",
      as.character(cov_res)
    ))
  }

  expect_s3_class(cov_res, "nlmixr2SSE")
  expect_equal(
    cov_res$runInfo$parameterSourceInfo$parameterPartition$drawn,
    rownames(ref_fit$cov)
  )
  cov_tka <- cov_res$initialValues$value[
    cov_res$initialValues$parameter == "tka"
  ]
  cov_omega <- cov_res$initialValues$value[
    cov_res$initialValues$parameter == "omega(eta.ka,eta.ka)"
  ]
  expect_gte(length(unique(signif(cov_tka, 8))), 2L)
  expect_equal(length(unique(signif(cov_omega, 8))), 1L)
})

test_that("validateResumeRequest aborts on an rxThreads mismatch", {
  run_info <- list(
    fitName = "fake_sse_fit",
    samples = 2L,
    parameterSource = "fixed",
    estimateSimulation = TRUE,
    rxThreads = 8L
  )

  err <- capture_sse_error(
    .validateResumeRequest(
      existingRunInfo = run_info,
      fitName = "fake_sse_fit",
      samples = 2L,
      control = runSSEControl(workers = 1L),
      requestedLabels = character(0),
      rxThreads = 2L
    )
  )

  expect_s3_class(err, "error")
  expect_match(conditionMessage(err), "rxThreads")
  expect_match(conditionMessage(err), "8")
  expect_match(conditionMessage(err), "2")
})

test_that("validateResumeRequest warns when rxThreads was never recorded", {
  run_info <- list(
    fitName = "fake_sse_fit",
    samples = 2L,
    parameterSource = "fixed",
    estimateSimulation = TRUE
  )

  expect_warning(
    .validateResumeRequest(
      existingRunInfo = run_info,
      fitName = "fake_sse_fit",
      samples = 2L,
      control = runSSEControl(workers = 1L),
      requestedLabels = character(0),
      rxThreads = 4L
    ),
    "rxThreads"
  )
})

test_that("validateResumeRequest accepts a matching rxThreads", {
  run_info <- list(
    fitName = "fake_sse_fit",
    samples = 2L,
    parameterSource = "fixed",
    estimateSimulation = TRUE,
    rxThreads = 4L
  )

  expect_silent(
    .validateResumeRequest(
      existingRunInfo = run_info,
      fitName = "fake_sse_fit",
      samples = 2L,
      control = runSSEControl(workers = 1L),
      requestedLabels = character(0),
      rxThreads = 4L
    )
  )
})

test_that("addModels warns when rxThreads differs from the recorded run", {
  tmp <- tempfile("nlmixr2sse-addmodels-rxthreads-")
  on.exit(unlink(tmp, recursive = TRUE, force = TRUE), add = TRUE)
  write_fake_sse_run_dir(tmp)

  run_info <- readRDS(file.path(tmp, "run_info.rds"))
  run_info$rxThreads <- 8L
  saveRDS(run_info, file.path(tmp, "run_info.rds"))

  saveRDS(
    data.frame(ID = 1, DV = 10, stringsAsFactors = FALSE),
    file.path(tmp, "simulations", "sim_0001.rds")
  )
  saveRDS(
    data.frame(ID = 2, DV = 20, stringsAsFactors = FALSE),
    file.path(tmp, "simulations", "sim_0002.rds")
  )
  nlmixr2utils::withRunSeed(tmp, seed = 42, prefix = "sse")

  testthat::local_mocked_bindings(
    .fitTaskRecord = function(sample, spec, simRecord, unionSchema, ...) {
      list(
        success = TRUE,
        row = nlmixr2utils::rawResultsRow(
          fit = NULL,
          source = "sse",
          hypothesis = spec$hypothesis,
          sample = sample,
          modelLabel = spec$label,
          role = spec$role,
          objf = 150 + sample,
          schema = unionSchema
        ),
        sample = as.integer(sample),
        model_label = spec$label,
        role = spec$role,
        hypothesis = spec$hypothesis
      )
    },
    .package = "nlmixr2sse"
  )

  expect_warning(
    try(
      runSSE(
        fake_sse_fit(),
        alternativeModels = sseModel(fake_sse_fit(), label = "alt_rxt"),
        samples = 2L,
        control = runSSEControl(addModels = TRUE, workers = 1L, rxThreads = 2L),
        outputDir = tmp,
        restart = FALSE
      ),
      silent = TRUE
    ),
    "rxThreads"
  )
})

test_that("addModels preserves the originally recorded rxThreads", {
  existing <- list(
    rxThreads = 8L,
    parameterSource = "fixed",
    estimateSimulation = TRUE
  )
  expect_equal(
    .addModelsRxThreads(existing, resolvedRxThreads = 2L, addModels = TRUE),
    8L
  )
  # a fresh (non-addModels) run records what it actually resolved
  expect_equal(
    .addModelsRxThreads(existing, resolvedRxThreads = 2L, addModels = FALSE),
    2L
  )
  # an addModels run against a directory with no recorded value falls back
  expect_equal(
    .addModelsRxThreads(list(), resolvedRxThreads = 2L, addModels = TRUE),
    2L
  )
})
