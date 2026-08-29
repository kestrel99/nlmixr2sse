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

test_that("joint draw incorporates the theta-omega covariance", {
  p <- 2L
  om0 <- matrix(c(0.30, 0.05, 0.05, 0.12), p, p)
  th0 <- c(tka = 0.45, add = 0.70)
  w0 <- .omegaToVec(om0)

  corTarget <- diag(5L)
  corTarget[1L, 3L] <- corTarget[3L, 1L] <- 0.45   # tka <-> om(1,1)
  corTarget[2L, 5L] <- corTarget[5L, 2L] <- -0.35  # add <-> om(2,2)
  sds <- c(0.18, 0.12, 0.060, 0.020, 0.030)
  sigma <- diag(sds) %*% corTarget %*% diag(sds)

  spec <- .jointDrawSpec(th0, list(list(omega = om0, index = seq_len(p))), sigma)

  set.seed(2024)
  draws <- lapply(seq_len(6000L), function(i) .drawJoint(spec))

  tka <- vapply(draws, function(d) unname(d$theta[["tka"]]), numeric(1))
  add <- vapply(draws, function(d) unname(d$theta[["add"]]), numeric(1))
  om11 <- vapply(draws, function(d) d$omega[[1L]][1L, 1L], numeric(1))
  om22 <- vapply(draws, function(d) d$omega[[1L]][2L, 2L], numeric(1))

  expect_equal(stats::cor(tka, om11), 0.45, tolerance = 0.05)
  expect_equal(stats::cor(add, om22), -0.35, tolerance = 0.05)
  # NOT the raw target SE of 0.060. om11 = L11^2 and phi1 = log(L11), so om11
  # is exactly lognormal with log-sd s = SE/omega = 0.2. The back-transform
  # therefore inflates the SD by a known, exact factor:
  #   E[om11]  = omega * exp(s^2 / 2)
  #   SD[om11] = E[om11] * sqrt(exp(s^2) - 1)  =  0.061829   (+3.05%)
  # Asserting 0.060 here is unpassable for ANY correct implementation. Assert
  # the closed-form value instead -- that keeps the test tight enough to catch
  # a real regression, where simply widening the tolerance would not.
  s <- 0.060 / 0.30
  sdClosedForm <- 0.30 * exp(s^2 / 2) * sqrt(exp(s^2) - 1)
  expect_equal(stats::sd(om11), sdClosedForm, tolerance = 0.03)
})

test_that("joint draw survives an OMEGA-only covariance (zero thetas)", {
  # regression guard for the -seq_len(0) trap: with no thetas, a negative-index
  # split would silently return an empty phi and lose the OMEGA draw entirely
  om0 <- matrix(c(0.30, 0.05, 0.05, 0.12), 2L, 2L)
  sigma <- diag(c(0.06, 0.02, 0.03)^2)

  spec <- .jointDrawSpec(
    stats::setNames(numeric(0), character(0)),
    list(list(omega = om0, index = 1:2)),
    sigma
  )

  set.seed(11)
  d <- .drawJoint(spec)

  expect_length(d$theta, 0L)
  expect_length(d$omega, 1L)
  expect_equal(dim(d$omega[[1L]]), c(2L, 2L))
  expect_true(all(is.finite(d$omega[[1L]])))
  expect_true(
    all(eigen(d$omega[[1L]], symmetric = TRUE, only.values = TRUE)$values > 0)
  )
})

test_that("every joint draw is positive-definite", {
  om0 <- matrix(c(0.30, 0.05, 0.05, 0.12), 2L, 2L)
  th0 <- c(tka = 0.45)
  sigma <- diag(c(0.03, 0.06, 0.02, 0.03)^2)

  spec <- .jointDrawSpec(th0, list(list(omega = om0, index = 1:2)), sigma)

  set.seed(5)
  ok <- vapply(seq_len(3000L), function(i) {
    om <- .drawJoint(spec)$omega[[1L]]
    all(eigen(om, symmetric = TRUE, only.values = TRUE)$values > 0)
  }, logical(1))
  expect_true(all(ok))
})

test_that("log-Cholesky round-trips across a wide scale range", {
  # OMEGA components span many orders of magnitude in practice. A fixed
  # absolute finite-difference step cannot serve all of them, which is why the
  # Jacobian step is relative.
  for (scale in c(1e-8, 1e-4, 1e-1, 1e1, 1e2)) {
    om <- matrix(c(1.0, 0.2, 0.2, 0.5), 2L, 2L) * scale
    phi <- .omegaToPhi(om)
    expect_equal(.phiToOmega(phi, 2L), om, tolerance = 1e-8,
                 info = paste("scale", scale))
  }
})

test_that("numericJacobian is full rank across a wide scale range", {
  for (scale in c(1e-8, 1e-4, 1e-1, 1e1, 1e2)) {
    om <- matrix(c(1.0, 0.2, 0.2, 0.5), 2L, 2L) * scale
    j <- .numericJacobian(
      function(w) .omegaToPhi(.vecToOmega(w, 2L)),
      .omegaToVec(om)
    )
    expect_true(all(is.finite(j)), info = paste("scale", scale))
    expect_equal(qr(j)$rank, 3L, info = paste("scale", scale))
  }
})

test_that("near-boundary correlations still transform and draw", {
  # correlation 0.98 -- valid but close to the edge of the PD cone
  om <- matrix(c(0.30, 0.98 * sqrt(0.30 * 0.12),
                 0.98 * sqrt(0.30 * 0.12), 0.12), 2L, 2L)
  expect_true(all(eigen(om, symmetric = TRUE, only.values = TRUE)$values > 0))

  phi <- .omegaToPhi(om)
  expect_equal(.phiToOmega(phi, 2L), om, tolerance = 1e-9)

  sigma <- diag(c(0.06, 0.02, 0.03)^2)
  spec <- .jointDrawSpec(
    stats::setNames(numeric(0), character(0)),
    list(list(omega = om, index = 1:2)),
    sigma
  )
  set.seed(21)
  ok <- vapply(seq_len(500L), function(i) {
    d <- .drawJoint(spec)$omega[[1L]]
    all(eigen(d, symmetric = TRUE, only.values = TRUE)$values > 0)
  }, logical(1))
  expect_true(all(ok))
})

test_that("an ill-conditioned but positive-definite block is handled", {
  om <- diag(c(1e-6, 1e2))   # condition number ~1e8
  expect_true(all(eigen(om, symmetric = TRUE, only.values = TRUE)$values > 0))
  j <- .numericJacobian(
    function(w) .omegaToPhi(.vecToOmega(w, 2L)),
    .omegaToVec(om)
  )
  expect_true(all(is.finite(j)))
  expect_equal(qr(j)$rank, 3L)
})

test_that("multi-block draws map back to the right etas", {
  # Two blocks of different sizes stacked into one joint draw. If the phi
  # vector were split at the wrong offsets, the blocks would be silently
  # transposed -- values would look plausible but belong to the wrong etas.
  omA <- matrix(c(0.30, 0.05, 0.05, 0.12), 2L, 2L)   # 3 phi elements
  omB <- matrix(9.00, 1L, 1L)                        # 1 phi element, distinct scale

  sigma <- diag(c(0.06, 0.02, 0.03, 1.50)^2)
  spec <- .jointDrawSpec(
    stats::setNames(numeric(0), character(0)),
    list(
      list(omega = omA, index = 1:2),
      list(omega = omB, index = 3L)
    ),
    sigma
  )

  set.seed(3)
  d <- .drawJoint(spec)

  expect_length(d$omega, 2L)
  expect_equal(dim(d$omega[[1L]]), c(2L, 2L))
  expect_equal(dim(d$omega[[2L]]), c(1L, 1L))
  # the second block's scale (~9) must not leak into the first (~0.3)
  expect_lt(d$omega[[1L]][1L, 1L], 3)
  expect_gt(d$omega[[2L]][1L, 1L], 3)
})
