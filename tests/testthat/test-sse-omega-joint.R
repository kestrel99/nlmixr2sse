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

test_that("numericJacobian matches a known analytic derivative", {
  # f(x) = (x1^2, 3*x2)  =>  J = [[2*x1, 0], [0, 3]]
  f <- function(x) c(x[[1L]]^2, 3 * x[[2L]])
  j <- .numericJacobian(f, c(2, 5))
  expect_equal(j, matrix(c(4, 0, 0, 3), 2L, 2L), tolerance = 1e-6)
})

test_that("numericJacobian handles the log-Cholesky map", {
  om0 <- matrix(c(0.30, 0.05, 0.05, 0.12), 2L, 2L)
  w0 <- .omegaToVec(om0)
  j <- .numericJacobian(function(w) .omegaToPhi(.vecToOmega(w, 2L)), w0)

  expect_equal(dim(j), c(3L, 3L))
  expect_true(all(is.finite(j)))
  # the map is invertible at a PD point, so the Jacobian must be full rank
  expect_equal(qr(j)$rank, 3L)
})

test_that("numericJacobian prefers a halved CENTRAL step over one-sided", {
  # f is undefined outside a narrow window at the base step but fine at h/2,
  # so pass 1's halving must succeed and pass 2 must NOT fire. A regression
  # that accepted one-sided at the first central failure would return a
  # first-order value here instead of the exact second-order one.
  x0 <- 1
  window <- 1e-6
  f <- function(x) {
    if (abs(x - x0) > window) stop("outside domain")
    x^2
  }
  j <- .numericJacobian(f, x0, typical = window * 4)
  expect_equal(as.numeric(j), 2 * x0, tolerance = 1e-6)
})

test_that("numericJacobian falls back to one-sided only when central is impossible", {
  # One side is permanently blocked at every step, so central can never
  # succeed and the one-sided fallback is the correct answer.
  f <- function(x) {
    if (x > 1) stop("upper side blocked")
    x^2
  }
  j <- .numericJacobian(f, 1, typical = 1e-3)
  expect_true(is.finite(as.numeric(j)))
  # one-sided is first-order, so a looser tolerance than the central case
  expect_equal(as.numeric(j), 2, tolerance = 1e-2)
})

test_that("numericJacobian aborts when no step works, naming the context", {
  f <- function(x) {
    if (length(x) && identical(x, 1)) return(1)
    stop("every perturbation fails")
  }
  err <- capture_sse_error(.numericJacobian(f, 1, context = "eta.a, eta.b"))
  expect_s3_class(err, "error")
  expect_match(conditionMessage(err), "eta\\.a")
})

test_that("phiToOmega rejects a wrong-length phi vector", {
  err <- capture_sse_error(.phiToOmega(c(1, 2), 2L))
  expect_s3_class(err, "error")
  expect_match(conditionMessage(err), "3")
})

test_that("vecToOmega rejects a wrong-length vector", {
  err <- capture_sse_error(.vecToOmega(c(1, 2), 2L))
  expect_s3_class(err, "error")
  expect_match(conditionMessage(err), "3")
})
