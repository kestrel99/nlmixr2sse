skip_on_cran()

test_that("runSSEControl returns a validated control object", {
  ctl <- runSSEControl()

  expect_s3_class(ctl, "nlmixr2SSEControl")
  expect_equal(ctl$parameterSource, "fixed")
  expect_equal(ctl$offsetRawres, 1L)
  expect_true(ctl$estimateSimulation)
})

test_that("runSSEControl enforces rawres-only options", {
  # The deprecated randomEstimationInits keeps its original hard restriction
  # to "rawres" (emitting a deprecation warning first); referenceInitials /
  # alternativeInitials, set directly, deliberately relax this to a no-op --
  # see the "simulation-start initials" tests below.
  err_random <- capture_sse_error(
    suppressWarnings(
      runSSEControl(
        parameterSource = "fixed",
        randomEstimationInits = TRUE
      )
    )
  )
  expect_s3_class(err_random, "error")
  expect_match(conditionMessage(err_random), "randomEstimationInits")

  err_update <- capture_sse_error(
    runSSEControl(parameterSource = "fixed", updateFix = TRUE)
  )
  expect_s3_class(err_update, "error")
  expect_match(conditionMessage(err_update), "updateFix")

  err_filter <- capture_sse_error(
    runSSEControl(
      parameterSource = "fixed",
      inFilter = "minimization_successful.eq.1"
    )
  )
  expect_s3_class(err_filter, "error")
  expect_match(conditionMessage(err_filter), "inFilter")
})

test_that("runSSEControl enforces staged-availability and combination rules", {
  err_ref <- capture_sse_error(
    runSSEControl(refOfv = 100, estimateSimulation = TRUE)
  )
  expect_s3_class(err_ref, "error")
  expect_match(conditionMessage(err_ref), "refOfv")

  err_append <- capture_sse_error(
    runSSEControl(appendColumns = "IPRED", addModels = TRUE)
  )
  expect_s3_class(err_append, "error")
  expect_match(conditionMessage(err_append), "appendColumns")

  ctl <- runSSEControl(simulationPostProcess = identity)
  expect_s3_class(ctl, "nlmixr2SSEControl")
  expect_identical(ctl$simulationPostProcess, identity)

  err_eta <- capture_sse_error(runSSEControl(initialEtas = TRUE))
  expect_s3_class(err_eta, "error")
  expect_match(conditionMessage(err_eta), "reserved")
})

test_that("initials are settable per model role", {
  # parameterSource = "rawres" avoids the "simulation has no effect" warning
  # covered separately below -- this test is only about the values sticking.
  ctl <- runSSEControl(
    parameterSource = "rawres",
    referenceInitials = "simulation",
    alternativeInitials = "model"
  )

  expect_equal(ctl$referenceInitials, "simulation")
  expect_equal(ctl$alternativeInitials, "model")
})

test_that("initials reject unknown policies", {
  expect_error(runSSEControl(referenceInitials = "truth"), "must be one of")
})

test_that("randomEstimationInits maps to both roles with a deprecation warning", {
  # The mapped policy is "simulation" for both roles, which -- like the
  # deprecated flag it replaces -- requires parameterSource = "rawres"; the
  # plan's sketch of this test omitted that and would abort on the default
  # "fixed" source instead of just warning, so it is supplied here to isolate
  # the deprecation-mapping behavior from that separate restriction.
  expect_warning(
    ctl <- runSSEControl(
      parameterSource = "rawres",
      randomEstimationInits = TRUE
    ),
    "deprecated"
  )
  expect_equal(ctl$referenceInitials, "simulation")
  expect_equal(ctl$alternativeInitials, "simulation")
})

test_that("the old and new arguments cannot contradict each other", {
  expect_error(
    suppressWarnings(runSSEControl(randomEstimationInits = TRUE,
                                   referenceInitials = "model")),
    "contradict"
  )
})

test_that("simulation-start initials warn (not error) outside rawres mode", {
  # Unlike the deprecated randomEstimationInits (a hard error),
  # referenceInitials / alternativeInitials = "simulation" outside "rawres"
  # warns: setting one role while the other keeps its default "model" is
  # legitimate, so it cannot be a hard error -- but every other mode-gated
  # argument in this function (updateFix, inFilter, covarianceDraw,
  # omegaRseWarn) never silently swallows an explicitly-set-but-inapplicable
  # option, so it cannot be silent either.
  # Each argument warns independently -- setting both at once (as here)
  # legitimately fires two warnings, so each expect_warning() below sets only
  # one to "simulation" at a time and leaves the other at its silent "model"
  # default; suppressWarnings() then covers the combined case explicitly, to
  # confirm both warnings fire together without leaking an uncaught one into
  # the suite's WARN count.
  expect_warning(
    ctl_ref <- runSSEControl(
      parameterSource = "fixed",
      referenceInitials = "simulation"
    ),
    "no effect"
  )
  expect_equal(ctl_ref$referenceInitials, "simulation")

  expect_warning(
    ctl_alt <- runSSEControl(
      parameterSource = "fixed",
      alternativeInitials = "simulation"
    ),
    "no effect"
  )
  expect_equal(ctl_alt$alternativeInitials, "simulation")

  warnings_both <- testthat::capture_warnings(
    ctl_fixed <- runSSEControl(
      parameterSource = "fixed",
      referenceInitials = "simulation",
      alternativeInitials = "simulation"
    )
  )
  expect_equal(length(warnings_both), 2L)
  expect_match(warnings_both[[1L]], "referenceInitials")
  expect_match(warnings_both[[2L]], "alternativeInitials")
  expect_equal(ctl_fixed$referenceInitials, "simulation")
  expect_equal(ctl_fixed$alternativeInitials, "simulation")

  # Covariance mode gets the same warning, but for a different underlying
  # reason (per-replicate vectors differ there, but simulation starts are
  # not wired up for it).
  expect_warning(
    runSSEControl(
      parameterSource = "covariance",
      referenceInitials = "simulation"
    ),
    "no effect"
  )

  # A default "model" policy must never warn.
  expect_silent(runSSEControl(parameterSource = "fixed"))
  expect_silent(
    runSSEControl(parameterSource = "fixed", referenceInitials = "model")
  )
})

test_that("runSSEControl validates raw-results filters", {
  ctl <- runSSEControl(
    parameterSource = "rawres",
    inFilter = "minimization_successful.eq.1",
    outFilter = ~ significant_digits > 3
  )

  expect_s3_class(ctl, "nlmixr2SSEControl")
  expect_equal(ctl$parameterSource, "rawres")
})

test_that("runSSEControl exposes rxThreads with an auto default", {
  ctl <- runSSEControl()
  expect_equal(ctl$rxThreads, "auto")

  expect_equal(runSSEControl(rxThreads = 1L)$rxThreads, 1L)
  expect_equal(runSSEControl(rxThreads = 8L)$rxThreads, 8L)
  expect_null(runSSEControl(rxThreads = NULL)$rxThreads)
})

test_that("runSSEControl rejects invalid rxThreads", {
  for (bad in list(0L, -1L, 2.5, "many", c(1L, 2L), NA_integer_)) {
    label <- paste("rxThreads =", paste(deparse(bad), collapse = " "))
    err <- capture_sse_error(runSSEControl(rxThreads = bad))
    expect_s3_class(err, "error")
    expect_match(conditionMessage(err), "rxThreads", info = label)
  }
})

test_that("runSSEControl exposes covarianceDraw defaulting to independent_iw", {
  expect_equal(runSSEControl()$covarianceDraw, "independent_iw")
  expect_equal(
    runSSEControl(
      parameterSource = "covariance",
      covarianceDraw = "joint"
    )$covarianceDraw,
    "joint"
  )
})

test_that("runSSEControl exposes omegaRseWarn", {
  expect_equal(runSSEControl()$omegaRseWarn, 0.5)
  expect_equal(
    runSSEControl(
      parameterSource = "covariance",
      omegaRseWarn = 0.25
    )$omegaRseWarn,
    0.25
  )
  err <- capture_sse_error(
    runSSEControl(parameterSource = "covariance", omegaRseWarn = -1)
  )
  expect_s3_class(err, "error")
})

test_that("runSSEControl rejects an unknown covarianceDraw", {
  err <- capture_sse_error(
    runSSEControl(parameterSource = "covariance", covarianceDraw = "wishart")
  )
  expect_s3_class(err, "error")
})

test_that("covarianceDraw requires parameterSource = covariance", {
  err <- capture_sse_error(
    runSSEControl(parameterSource = "fixed", covarianceDraw = "independent_iw")
  )
  expect_s3_class(err, "error")
  expect_match(conditionMessage(err), "covarianceDraw")
})

test_that("workerDescription reports the resolved rxode2 thread count", {
  expect_match(
    nlmixr2sse:::.workerDescription(1L, rxThreads = 8L),
    "sequential \\(workers = 1\\), 8 rxode2 thread"
  )
  expect_match(
    nlmixr2sse:::.workerDescription(4L, rxThreads = 2L),
    "parallel \\(4 workers\\), 2 rxode2 thread"
  )
  # omitted thread count keeps the bare description
  expect_equal(
    nlmixr2sse:::.workerDescription(1L),
    "sequential (workers = 1)"
  )
})
