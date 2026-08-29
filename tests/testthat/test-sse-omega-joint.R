skip_on_cran()

test_that("log-Cholesky transform round-trips", {
  for (om in list(
    matrix(0.3, 1L, 1L),
    matrix(c(0.30, 0.05, 0.05, 0.12), 2L, 2L),
    matrix(c(0.40, 0.05, 0.02,
             0.05, 0.20, 0.03,
             0.02, 0.03, 0.15), 3L, 3L)
  )) {
    phi <- .omegaToPhi(om)
    expect_equal(.phiToOmega(phi, nrow(om)), om, tolerance = 1e-10)
  }
})

test_that("phiToOmega is positive-definite for arbitrary input", {
  set.seed(3)
  for (i in seq_len(200L)) {
    phi <- stats::rnorm(3L, sd = 5)          # deliberately wild values
    om <- .phiToOmega(phi, 2L)
    expect_true(all(eigen(om, symmetric = TRUE, only.values = TRUE)$values > 0))
  }
})

test_that("omegaToPhi logs the diagonal only", {
  om <- diag(c(exp(2), exp(4)))
  # L = diag(exp(1), exp(2)); phi = (log L11, L21, log L22) = (1, 0, 2)
  expect_equal(.omegaToPhi(om), c(1, 0, 2))
})
