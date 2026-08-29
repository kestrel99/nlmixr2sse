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
#' positive-definite for any `phi`.
#' @noRd
.phiToOmega <- function(phi, p) {
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
  idx <- .lowerTriIndex(p)
  omega <- matrix(0, p, p)
  omega[cbind(idx[, "row"], idx[, "col"])] <- v
  omega[cbind(idx[, "col"], idx[, "row"])] <- v
  omega
}
