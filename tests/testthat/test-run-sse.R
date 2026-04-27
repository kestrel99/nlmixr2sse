skip_on_cran()

test_that("runSSE initializes a numbered run directory and returns the SSE class", {
  tmp <- tempfile("nlmixr2sse-phase3-")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE, force = TRUE), add = TRUE)
  old <- setwd(tmp)
  on.exit(setwd(old), add = TRUE)

  fit <- fake_sse_fit()
  res <- runSSE(fit, samples = 5L, restart = TRUE)

  expect_s3_class(res, "nlmixr2SSE")
  expect_equal(basename(res$outputDir), "fit_sse_1")
  expect_true(file.exists(file.path(res$outputDir, "run_info.rds")))
  expect_true(file.exists(file.path(res$outputDir, "sse_seed.rds")))
  expect_true(dir.exists(file.path(res$outputDir, "simulations")))
  expect_true(dir.exists(file.path(res$outputDir, "fits")))
  expect_equal(nrow(res$rawResults), 0L)
  expect_equal(res$runInfo$status, "initialized")
})

test_that("runSSE derives fitName from expressions and resolves alternative labels", {
  tmp <- tempfile("nlmixr2sse-explicit-")
  on.exit(unlink(tmp, recursive = TRUE, force = TRUE), add = TRUE)

  fits <- list(ref = fake_sse_fit())
  alt1 <- sseModel(function() NULL, est = "focei", control = list(print = 0L))
  alt2 <- sseModel(
    function() NULL,
    est = "saem",
    control = list(print = 0L),
    label = "reduced"
  )

  res <- runSSE(
    fits$ref,
    alternativeModels = list(alt1, alt2),
    samples = 2L,
    outputDir = tmp,
    restart = TRUE
  )

  expect_equal(res$runInfo$fitName, "fits_ref")
  expect_equal(
    vapply(res$alternativeSpecs, `[[`, character(1), "label"),
    c("alt1", "reduced")
  )
})

test_that("runSSE validates covariance mode and alternative-label collisions", {
  bad_fit <- fake_sse_fit()
  bad_fit$cov <- NULL

  err_cov <- capture_sse_error(
    runSSE(
      bad_fit,
      control = runSSEControl(parameterSource = "covariance"),
      outputDir = tempfile("nlmixr2sse-badcov-"),
      restart = TRUE
    )
  )
  expect_s3_class(err_cov, "error")
  expect_match(conditionMessage(err_cov), "fit\\$cov")

  err_label <- capture_sse_error(
    runSSE(
      fake_sse_fit(),
      alternativeModels = sseModel(
        function() NULL,
        est = "focei",
        control = list(print = 0L),
        label = "fake_sse_fit"
      ),
      outputDir = tempfile("nlmixr2sse-badlabel-"),
      restart = TRUE
    )
  )
  expect_s3_class(err_label, "error")
  expect_match(conditionMessage(err_label), "collides")
})
