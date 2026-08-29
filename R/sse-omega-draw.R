# Inverse-Wishart machinery for drawing OMEGA with uncertainty.
#
# Shares NWPRI's broad independence factorization -- THETA and OMEGA drawn from
# INDEPENDENT distributions (Normal and inverse-Wishart respectively), so the
# THETA<->OMEGA cross-terms present in fit$cov are deliberately not used here.
# This is NOT NONMEM $PRIOR NWPRI: the OMEGA density differs.

#' Name one OMEGA entry the way nlmixr2est names it in `fit$cov`
#'
#' Diagonal entries are `om.<eta>`; off-diagonal entries are
#' `cov.<rowEta>.<colEta>`.
#' @noRd
.omegaCovName <- function(rowName, colName) {
  ifelse(
    rowName == colName,
    paste0("om.", rowName),
    paste0("cov.", rowName, ".", colName)
  )
}

#' Tabulate OMEGA entries with their `fit$cov` names
#'
#' @param info a `.uiOmegaInfo()` result (columns row, col, rowName, colName,
#'   fix). Remember `.uiOmegaInfo()` must be given `fit$ui`, not `fit`.
#' @return data frame with row, col, rowName, colName, fix, covName, diagonal
#' @noRd
.omegaEntryTable <- function(info) {
  if (!is.data.frame(info) || nrow(info) == 0L) {
    return(data.frame(
      row = integer(0),
      col = integer(0),
      rowName = character(0),
      colName = character(0),
      fix = logical(0),
      covName = character(0),
      diagonal = logical(0),
      stringsAsFactors = FALSE
    ))
  }

  data.frame(
    row = as.integer(info$row),
    col = as.integer(info$col),
    rowName = as.character(info$rowName),
    colName = as.character(info$colName),
    fix = as.logical(info$fix),
    covName = .omegaCovName(
      as.character(info$rowName),
      as.character(info$colName)
    ),
    diagonal = as.integer(info$row) == as.integer(info$col),
    stringsAsFactors = FALSE
  )
}
