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

#' Split OMEGA into independent blocks
#'
#' Two etas belong to the same block when an off-diagonal entry links them,
#' directly or transitively. A fully diagonal OMEGA yields all 1x1 blocks.
#'
#' @param tab an `.omegaEntryTable()` result
#' @param nEta total number of etas
#' @return list of integer vectors, each sorted, ordered by first member
#' @noRd
.omegaBlocks <- function(tab, nEta) {
  nEta <- as.integer(nEta)
  if (nEta <= 0L) {
    return(list())
  }

  # union-find over eta indices
  parent <- seq_len(nEta)
  findRoot <- function(i) {
    while (parent[[i]] != i) {
      i <- parent[[i]]
    }
    i
  }
  offDiag <- tab[!tab$diagonal, , drop = FALSE]
  for (i in seq_len(nrow(offDiag))) {
    a <- findRoot(as.integer(offDiag$row[[i]]))
    b <- findRoot(as.integer(offDiag$col[[i]]))
    if (a != b) {
      parent[[max(a, b)]] <- min(a, b)
    }
  }

  roots <- vapply(seq_len(nEta), findRoot, integer(1))
  # order blocks by their first member so output is deterministic
  unname(split(seq_len(nEta), factor(roots, levels = unique(roots))))
}

#' Apply the OMEGA coverage policy
#'
#' A block is drawable only when EVERY declared entry in it is unfixed and
#' present in `fit$cov`. A block with any fixed or uncovered element is held
#' entirely at its fitted values.
#'
#' This is deliberately stricter than "some diagonal has a usable standard
#' error". `nlmixr2est` drops fixed OMEGA elements from `fit$cov` while the
#' declared topology still lists them, so a correlated block with one fixed
#' component has missing entries. Drawing such a block would mutate the fixed
#' element; drawing the free sub-matrix and splicing the fixed values back in
#' can produce a non-positive-definite matrix. Neither is acceptable, so the
#' whole block is held.
#'
#' Both draw modes MUST consume this result rather than re-deriving
#' drawability, so that `"joint"` and `"independent_iw"` agree on what varies.
#'
#' @param blocks list of integer eta-index vectors, from `.omegaBlocks()`
#' @param entries an `.omegaEntryTable()` result
#' @param covNames `rownames(fit$cov)`
#' @return list(drawable = <list of index vectors>,
#'              held = <list of list(index, reason)>)
#' @noRd
.drawableOmegaBlocks <- function(blocks, entries, covNames) {
  drawable <- list()
  held <- list()

  for (idx in blocks) {
    inBlock <- entries$row %in% idx & entries$col %in% idx
    blockEntries <- entries[inBlock, , drop = FALSE]

    if (nrow(blockEntries) == 0L) {
      .abortSSE(
        "Internal error: OMEGA block {.val {paste(idx, collapse = ', ')}} has no entries in the entry table."
      )
    }

    fixedNames <- blockEntries$covName[blockEntries$fix]
    missingNames <- setdiff(
      blockEntries$covName[!blockEntries$fix],
      covNames
    )

    if (length(fixedNames) == 0L && length(missingNames) == 0L) {
      drawable[[length(drawable) + 1L]] <- idx
      next
    }

    reason <- paste(
      c(
        if (length(fixedNames) > 0L) {
          paste0("fixed: ", paste(fixedNames, collapse = ", "))
        },
        if (length(missingNames) > 0L) {
          paste0("not covered by fit$cov: ", paste(missingNames, collapse = ", "))
        }
      ),
      collapse = "; "
    )
    held[[length(held) + 1L]] <- list(index = idx, reason = reason)
  }

  list(drawable = drawable, held = held)
}
