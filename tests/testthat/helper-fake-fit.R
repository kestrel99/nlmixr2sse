`$.nlmixr2SSEFakeFit` <- function(x, name) {
  x[[name]]
}

fake_sse_fit <- function() {
  fit <- list(
    finalUiEnv = function() NULL,
    est = "focei",
    control = list(print = 0L),
    theta = c(tka = 0.45, tcl = 1),
    omega = matrix(
      0.3,
      1L,
      1L,
      dimnames = list("eta.ka", "eta.ka")
    ),
    sigma = matrix(numeric(0), 0L, 0L),
    parFixedDf = data.frame(
      Estimate = c(0.45, 1),
      `Back-transformed` = c(exp(0.45), exp(1)),
      SE = c(0.1, 0.2),
      row.names = c("tka", "tcl"),
      check.names = FALSE
    ),
    iniDf = data.frame(
      name = c("tka", "tcl", "eta.ka"),
      ntheta = c(1L, 2L, NA),
      neta1 = c(NA, NA, 1L),
      neta2 = c(NA, NA, 1L),
      lower = c(-Inf, -Inf, NA),
      upper = c(Inf, Inf, NA),
      fix = c(FALSE, FALSE, FALSE),
      stringsAsFactors = FALSE
    ),
    cov = matrix(
      c(0.04, 0.01, 0.01, 0.09),
      2L,
      2L,
      dimnames = list(c("tka", "tcl"), c("tka", "tcl"))
    ),
    objf = 123.45
  )
  class(fit) <- c("nlmixr2SSEFakeFit", "nlmixr2FitCore")
  fit
}

capture_sse_error <- function(expr) {
  tryCatch(
    {
      force(expr)
      NULL
    },
    error = function(e) e
  )
}
