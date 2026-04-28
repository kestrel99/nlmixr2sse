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
