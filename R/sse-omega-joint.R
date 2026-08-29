# Joint THETA+OMEGA draw on a log-Cholesky-transformed scale.
#
# Unlike the "independent_iw" mode, this incorporates the THETA<->OMEGA covariance that
# fit$cov reports. OMEGA is transformed to an unconstrained vector so that a
# plain multivariate-Normal draw always back-transforms to a positive-definite
# matrix -- no rejection sampling.

#' Lower-triangular index pairs, column-major, diagonal included
#' @noRd
.lowerTriIndex <- function(p) {
  idx <- which(lower.tri(matrix(0, p, p), diag = TRUE), arr.ind = TRUE)
  idx[order(idx[, "col"], idx[, "row"]), , drop = FALSE]
}

#' Map an OMEGA block to its unconstrained log-Cholesky vector
#' @noRd
.omegaToPhi <- function(omega) {
  p <- nrow(omega)
  lower <- t(chol(omega))
  idx <- .lowerTriIndex(p)
  v <- lower[cbind(idx[, "row"], idx[, "col"])]
  onDiag <- idx[, "row"] == idx[, "col"]
  v[onDiag] <- log(v[onDiag])
  v
}

#' Rebuild an OMEGA block from its log-Cholesky vector
#'
#' Exponentiating the diagonal keeps it strictly positive, so the result is
#' positive-definite by construction in exact arithmetic (L is triangular with
#' strictly positive diagonal). At extreme condition numbers -- around 1e16,
#' reachable with deliberately wild synthetic input on 3x3+ blocks --
#' floating-point roundoff can still produce a matrix that fails a strict
#' eigen() check. Realistic delta-method draws centred on a fitted OMEGA do
#' not reach that regime, but callers should not assume ironclad PD for large
#' ill-conditioned blocks.
#' @noRd
.phiToOmega <- function(phi, p) {
  expected <- p * (p + 1L) / 2L
  if (length(phi) != expected) {
    .abortSSE(
      "Expected {.val {expected}} log-Cholesky elements for a {.val {p}}-eta block, got {.val {length(phi)}}."
    )
  }
  idx <- .lowerTriIndex(p)
  onDiag <- idx[, "row"] == idx[, "col"]
  v <- phi
  v[onDiag] <- exp(v[onDiag])
  lower <- matrix(0, p, p)
  lower[cbind(idx[, "row"], idx[, "col"])] <- v
  lower %*% t(lower)
}

#' Lower-triangular elements of a symmetric matrix, matching `.lowerTriIndex()`
#' @noRd
.omegaToVec <- function(omega) {
  idx <- .lowerTriIndex(nrow(omega))
  omega[cbind(idx[, "row"], idx[, "col"])]
}

#' Rebuild a symmetric matrix from `.omegaToVec()` output
#' @noRd
.vecToOmega <- function(v, p) {
  expected <- p * (p + 1L) / 2L
  if (length(v) != expected) {
    .abortSSE(
      "Expected {.val {expected}} lower-triangular elements for a {.val {p}}-eta block, got {.val {length(v)}}."
    )
  }
  idx <- .lowerTriIndex(p)
  omega <- matrix(0, p, p)
  omega[cbind(idx[, "row"], idx[, "col"])] <- v
  omega[cbind(idx[, "col"], idx[, "row"])] <- v
  omega
}

#' Scale-aware central-difference Jacobian of a vector-valued function
#'
#' Used to carry `fit$cov` from the natural OMEGA scale onto the log-Cholesky
#' scale.
#'
#' The step MUST be relative, not absolute. OMEGA components span many orders
#' of magnitude (roughly 1e-8 to 1e2 in practice); a fixed `h = 1e-6` is too
#' coarse for small components, too fine for large ones, and can push a
#' near-boundary covariance matrix outside the positive-definite cone, where
#' the transform's `chol()` fails outright. So: relative step, halve on
#' failure, fall back to one-sided, and abort rather than return a silently
#' wrong derivative.
#' @noRd
.numericJacobian <- function(f, x, typical = 1e-4, maxHalve = 8L,
                             context = NULL) {
  f0 <- f(x)
  out <- matrix(0, length(f0), length(x))
  eps13 <- .Machine$double.eps^(1 / 3)

  for (j in seq_along(x)) {
    h0 <- max(abs(x[[j]]), typical) * eps13
    col <- NULL

    # PASS 1: exhaust central differencing, halving the step each time.
    # One-sided differences are first-order accurate where central is
    # second-order, so a one-sided result is a genuine loss of accuracy and
    # must not be accepted merely because the FIRST central step failed.
    h <- h0
    lastUp <- NULL
    lastDown <- NULL
    lastH <- h0
    for (attempt in seq_len(maxHalve)) {
      up <- x
      down <- x
      up[[j]] <- up[[j]] + h
      down[[j]] <- down[[j]] - h

      fUp <- tryCatch(f(up), error = function(e) NULL)
      fDown <- tryCatch(f(down), error = function(e) NULL)

      if (!is.null(fUp) && !is.null(fDown)) {
        col <- (fUp - fDown) / (2 * h)
        break
      }
      # remember the smallest step where at least one side worked, for pass 2
      if (!is.null(fUp) || !is.null(fDown)) {
        lastUp <- fUp
        lastDown <- fDown
        lastH <- h
      }
      h <- h / 2
    }

    # PASS 2: only now, having failed central differencing at every step,
    # fall back to one-sided on whichever side remained admissible.
    if (is.null(col) && !is.null(lastUp)) {
      col <- (lastUp - f0) / lastH
    }
    if (is.null(col) && !is.null(lastDown)) {
      col <- (f0 - lastDown) / lastH
    }

    if (is.null(col) || !all(is.finite(col))) {
      .abortSSE(
        paste0(
          "Could not differentiate the OMEGA transform",
          if (!is.null(context)) " for block {.val {context}}" else "",
          " at element {.val {j}}: every trial step left the ",
          "positive-definite cone. The fitted OMEGA is probably near-singular."
        )
      )
    }
    out[, j] <- col
  }
  out
}
