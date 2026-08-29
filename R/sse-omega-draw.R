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

#' Build the inverse-Wishart parameters for one OMEGA block
#'
#' Centres the distribution on `omega0` and matches the spread to the reported
#' OMEGA standard errors. For `Omega ~ InvWishart(Psi, nu)` over a p x p block,
#' setting `Psi = (nu - p - 1) * omega0` gives `E[Omega] = omega0` and
#' `Var(Omega_ii) = 2 * omega0_ii^2 / (nu - p - 3)`, so
#' `nu_i = p + 3 + 2 * (omega0_ii / se_i)^2`.
#'
#' The inverse-Wishart has a single `nu` per block, so only the binding
#' (minimum-nu) element matches its SE exactly; every other element comes out
#' MORE dispersed than reported. That is wider for the constructed diagonal
#' marginal, not a general conservatism guarantee, and is
#' intentional.
#'
#' @param omega0 the block's fitted OMEGA sub-matrix
#' @param se reported standard errors for the block's diagonal, `NA` where
#'   unavailable
#' @return list(omega0, nu, p), or `NULL` when no diagonal element has a usable SE
#' @noRd
.omegaWishartSpec <- function(omega0, se) {
  p <- nrow(omega0)
  variances <- diag(omega0)

  usable <- !is.na(se) & is.finite(se) & se > 0
  if (!any(usable)) {
    return(NULL)
  }

  if (any(!is.finite(variances[usable]) | variances[usable] <= 0)) {
    bad <- which(usable & (!is.finite(variances) | variances <= 0))
    badLabel <- if (!is.null(rownames(omega0))) {
      rownames(omega0)[bad]
    } else {
      as.character(bad)
    }
    .abortSSE(
      paste0(
        "OMEGA variance {.val {badLabel}} is zero or non-finite, so its uncertainty ",
        "cannot be characterised. Fix the parameter, or use ",
        "{.code parameterSource = \"fixed\"}."
      )
    )
  }

  nuEach <- p + 3 + 2 * (variances[usable] / se[usable])^2
  nu <- min(nuEach)

  # nu = p + 3 + (positive term) by construction, so `nu <= p + 3` catches only
  # underflow or invalid input -- NOT ordinary weak identification. Weak
  # identification is surfaced separately, by the relative-SE warning in
  # .warnWeakOmega(); do not try to detect it here.
  if (!is.finite(nu)) {
    .abortSSE(
      paste0(
        "The reported OMEGA standard errors imply a non-finite inverse-Wishart ",
        "degrees of freedom for a {.val {p}}-eta block."
      )
    )
  }

  ev <- suppressWarnings(
    eigen(omega0, symmetric = TRUE, only.values = TRUE)$values
  )
  if (!all(is.finite(ev)) || min(ev) <= 0) {
    .abortSSE(
      paste0(
        "A fitted OMEGA block is not positive-definite, so it cannot be used ",
        "as an inverse-Wishart centre."
      )
    )
  }

  list(omega0 = omega0, nu = nu, p = p)
}

#' Draw one OMEGA block from a mean-centred inverse-Wishart
#'
#' Delegates the draw to `rxode2::cvPost()` rather than re-implementing it.
#'
#' `cvPost(nu, omega)` uses the scale convention `Psi = nu * omega`, so it is
#' NOT mean-centred: `E[Omega*] = nu/(nu - p - 1) * omega`. Verified empirically
#' -- at `nu = 20, p = 2` the observed ratio is 1.179 against a predicted 1.176.
#' Pre-scaling the input by `(nu - p - 1)/nu` cancels that exactly, recovering
#' `E[Omega*] = omega` and the textbook `Var = 2*omega_ii^2/(nu - p - 3)` that
#' `.omegaWishartSpec()`'s moment match assumes.
#'
#' The result is positive-definite by construction, so no rejection sampling is
#' needed.
#' @noRd
.drawOmegaBlock <- function(spec) {
  scaled <- spec$omega0 * (spec$nu - spec$p - 1) / spec$nu
  drawn <- tryCatch(
    rxode2::cvPost(spec$nu, scaled, n = 1L),
    error = function(e) {
      .abortSSE(
        "An OMEGA block could not be drawn from its inverse-Wishart: {conditionMessage(e)}"
      )
    }
  )
  if (is.list(drawn)) {
    drawn <- drawn[[1L]]
  }
  # enforce exact symmetry (guards against accumulated numerical asymmetry)
  0.5 * (drawn + t(drawn))
}

#' Draw a full OMEGA matrix, block by block
#'
#' `blocks` MUST already be `.drawableOmegaBlocks()$drawable` — this function
#' does not decide drawability, and must not be handed the raw block list.
#' Anything not in `blocks` keeps its fitted value.
#'
#' @param omega0 the fitted OMEGA
#' @param blocks DRAWABLE eta-index vectors only
#' @param se per-eta reported standard errors
#' @return a full OMEGA matrix with `omega0`'s dimnames
#' @noRd
.drawOmega <- function(omega0, blocks, se) {
  out <- omega0
  for (idx in blocks) {
    spec <- .omegaWishartSpec(
      omega0[idx, idx, drop = FALSE],
      se[idx]
    )
    if (is.null(spec)) {
      # a drawable block with no usable SE cannot be given a nu; keep it fitted
      next
    }
    out[idx, idx] <- .drawOmegaBlock(spec)
  }
  out
}
