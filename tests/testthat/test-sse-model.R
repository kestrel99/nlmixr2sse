skip_on_cran()

test_that("sseModel uses fitted-model defaults", {
  fit <- fake_sse_fit()

  spec <- sseModel(fit)

  expect_s3_class(spec, "nlmixr2SSEModel")
  expect_true(spec$isFit)
  expect_equal(spec$est, "focei")
  expect_equal(spec$control, fit[["control"]])
  expect_null(spec$label)
})

test_that("sseModel requires est and control for unfitted models", {
  mod <- function() NULL

  err_est <- capture_sse_error(sseModel(mod, control = list(print = 0L)))
  expect_s3_class(err_est, "error")
  expect_match(conditionMessage(err_est), "Supply")

  err_ctl <- capture_sse_error(sseModel(mod, est = "focei"))
  expect_s3_class(err_ctl, "error")
  expect_match(conditionMessage(err_ctl), "control")
})

test_that("sseModel validates labels", {
  spec <- sseModel(
    fake_sse_fit(),
    label = "altA"
  )

  expect_equal(spec$label, "altA")

  err <- capture_sse_error(
    sseModel(
      fake_sse_fit(),
      label = "bad label"
    )
  )
  expect_s3_class(err, "error")
  expect_match(conditionMessage(err), "letters")
})

test_that("mergeRawResultsSchemas emits no columns for a parameterless schema", {
  # paste0() drops zero-length arguments instead of returning zero length, so
  # paste0(character(0), ".se") is ".se". Without a guard, merging schemas that
  # contribute no parameters produced a spurious ".se" column, and rows built
  # from that schema had a column count no other schema could rbind against.
  empty <- list(
    thetaCols = character(0),
    omegaCols = character(0),
    sigmaCols = character(0)
  )

  merged <- .mergeRawResultsSchemas(list(empty, empty))

  expect_equal(merged$seCols, character(0))
  expect_false(".se" %in% merged$columns)
  expect_equal(merged$columns, merged$baseCols)
})

test_that("mergeRawResultsSchemas still derives se columns when parameters exist", {
  merged <- .mergeRawResultsSchemas(list(
    list(
      thetaCols = c("tka", "tcl"),
      omegaCols = "omega(eta.ka,eta.ka)",
      sigmaCols = character(0)
    )
  ))

  expect_equal(
    merged$seCols,
    c("tka.se", "tcl.se", "omega(eta.ka,eta.ka).se")
  )
  expect_true(all(merged$seCols %in% merged$columns))
})
