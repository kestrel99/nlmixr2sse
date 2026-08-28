skip_on_cran()

test_that("runSSEControl returns a validated control object", {
  ctl <- runSSEControl()

  expect_s3_class(ctl, "nlmixr2SSEControl")
  expect_equal(ctl$parameterSource, "fixed")
  expect_equal(ctl$offsetRawres, 1L)
  expect_true(ctl$estimateSimulation)
})

test_that("runSSEControl enforces rawres-only options", {
  err_random <- capture_sse_error(
    runSSEControl(
      parameterSource = "fixed",
      randomEstimationInits = TRUE
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
