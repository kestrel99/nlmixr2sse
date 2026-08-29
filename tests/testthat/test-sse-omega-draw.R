skip_on_cran()

test_that("omegaCovName builds nlmixr2est-style names", {
  expect_equal(.omegaCovName("eta.ka", "eta.ka"), "om.eta.ka")
  expect_equal(.omegaCovName("eta.cl", "eta.ka"), "cov.eta.cl.eta.ka")
})

test_that("omegaEntryTable maps positions to names and fixed flags", {
  info <- data.frame(
    row = c(1L, 2L, 2L, 3L),
    col = c(1L, 1L, 2L, 3L),
    rowName = c("eta.ka", "eta.cl", "eta.cl", "eta.v"),
    colName = c("eta.ka", "eta.ka", "eta.cl", "eta.v"),
    fix = c(FALSE, FALSE, FALSE, TRUE),
    stringsAsFactors = FALSE
  )

  tab <- .omegaEntryTable(info)

  expect_equal(
    tab$covName,
    c("om.eta.ka", "cov.eta.cl.eta.ka", "om.eta.cl", "om.eta.v")
  )
  expect_equal(tab$diagonal, c(TRUE, FALSE, TRUE, TRUE))
  expect_equal(tab$fix, c(FALSE, FALSE, FALSE, TRUE))
})

test_that("omegaEntryTable handles an empty omega", {
  info <- data.frame(
    row = integer(0), col = integer(0),
    rowName = character(0), colName = character(0),
    fix = logical(0), stringsAsFactors = FALSE
  )
  tab <- .omegaEntryTable(info)
  expect_equal(nrow(tab), 0L)
  expect_true(is.character(tab$covName))
})

test_that("omegaBlocks finds connected components", {
  # eta1+eta2 correlated; eta3 alone
  tab <- data.frame(
    row = c(1L, 2L, 2L, 3L),
    col = c(1L, 1L, 2L, 3L),
    diagonal = c(TRUE, FALSE, TRUE, TRUE),
    stringsAsFactors = FALSE
  )
  expect_equal(.omegaBlocks(tab, nEta = 3L), list(c(1L, 2L), 3L))
})

test_that("omegaBlocks treats a fully diagonal omega as 1x1 blocks", {
  tab <- data.frame(
    row = c(1L, 2L, 3L),
    col = c(1L, 2L, 3L),
    diagonal = c(TRUE, TRUE, TRUE),
    stringsAsFactors = FALSE
  )
  expect_equal(.omegaBlocks(tab, nEta = 3L), list(1L, 2L, 3L))
})

test_that("omegaBlocks merges transitively linked etas", {
  # 1-2 linked, 2-3 linked => one block of all three
  tab <- data.frame(
    row = c(1L, 2L, 2L, 3L, 3L),
    col = c(1L, 1L, 2L, 2L, 3L),
    diagonal = c(TRUE, FALSE, TRUE, FALSE, TRUE),
    stringsAsFactors = FALSE
  )
  expect_equal(.omegaBlocks(tab, nEta = 3L), list(c(1L, 2L, 3L)))
})

test_that("omegaBlocks handles zero etas", {
  tab <- data.frame(
    row = integer(0), col = integer(0), diagonal = logical(0),
    stringsAsFactors = FALSE
  )
  expect_equal(.omegaBlocks(tab, nEta = 0L), list())
})

test_that("drawableOmegaBlocks requires every declared entry to be covered", {
  entries <- data.frame(
    row = c(1L, 2L, 2L, 3L),
    col = c(1L, 1L, 2L, 3L),
    rowName = c("eta.a", "eta.b", "eta.b", "eta.c"),
    colName = c("eta.a", "eta.a", "eta.b", "eta.c"),
    fix = c(FALSE, FALSE, FALSE, FALSE),
    covName = c("om.eta.a", "cov.eta.b.eta.a", "om.eta.b", "om.eta.c"),
    diagonal = c(TRUE, FALSE, TRUE, TRUE),
    stringsAsFactors = FALSE
  )
  blocks <- list(c(1L, 2L), 3L)

  # everything covered -> both blocks drawable
  res <- .drawableOmegaBlocks(
    blocks, entries,
    covNames = c("om.eta.a", "cov.eta.b.eta.a", "om.eta.b", "om.eta.c")
  )
  expect_equal(res$drawable, list(c(1L, 2L), 3L))
  expect_length(res$held, 0L)

  # off-diagonal missing -> the WHOLE correlated block is held, not part of it
  res2 <- .drawableOmegaBlocks(
    blocks, entries,
    covNames = c("om.eta.a", "om.eta.b", "om.eta.c")
  )
  expect_equal(res2$drawable, list(3L))
  expect_equal(res2$held[[1L]]$index, c(1L, 2L))
  expect_match(res2$held[[1L]]$reason, "cov.eta.b.eta.a")
})

test_that("drawableOmegaBlocks holds a block containing a fixed element", {
  entries <- data.frame(
    row = c(1L, 2L, 2L),
    col = c(1L, 1L, 2L),
    rowName = c("eta.a", "eta.b", "eta.b"),
    colName = c("eta.a", "eta.a", "eta.b"),
    fix = c(FALSE, FALSE, TRUE),
    covName = c("om.eta.a", "cov.eta.b.eta.a", "om.eta.b"),
    diagonal = c(TRUE, FALSE, TRUE),
    stringsAsFactors = FALSE
  )
  res <- .drawableOmegaBlocks(
    list(c(1L, 2L)), entries,
    covNames = c("om.eta.a", "cov.eta.b.eta.a")
  )
  expect_length(res$drawable, 0L)
  expect_match(res$held[[1L]]$reason, "fixed")
})

test_that("drawableOmegaBlocks handles a theta-only covariance", {
  entries <- data.frame(
    row = 1L, col = 1L, rowName = "eta.a", colName = "eta.a",
    fix = FALSE, covName = "om.eta.a", diagonal = TRUE,
    stringsAsFactors = FALSE
  )
  res <- .drawableOmegaBlocks(list(1L), entries, covNames = c("tka", "tcl"))
  expect_length(res$drawable, 0L)
  expect_equal(res$held[[1L]]$index, 1L)
})

test_that("drawableOmegaBlocks reports both fixed and missing entries together", {
  entries <- data.frame(
    row = c(1L, 2L, 2L),
    col = c(1L, 1L, 2L),
    rowName = c("eta.a", "eta.b", "eta.b"),
    colName = c("eta.a", "eta.a", "eta.b"),
    fix = c(TRUE, FALSE, FALSE),
    covName = c("om.eta.a", "cov.eta.b.eta.a", "om.eta.b"),
    diagonal = c(TRUE, FALSE, TRUE),
    stringsAsFactors = FALSE
  )
  # om.eta.a is fixed AND cov.eta.b.eta.a is absent from covNames
  res <- .drawableOmegaBlocks(
    list(c(1L, 2L)), entries,
    covNames = "om.eta.b"
  )

  expect_length(res$drawable, 0L)
  expect_match(res$held[[1L]]$reason, "fixed: om\\.eta\\.a")
  expect_match(res$held[[1L]]$reason, "cov\\.eta\\.b\\.eta\\.a")
})

test_that("drawableOmegaBlocks aborts on a block/entries desync", {
  entries <- data.frame(
    row = 1L, col = 1L, rowName = "eta.a", colName = "eta.a",
    fix = FALSE, covName = "om.eta.a", diagonal = TRUE,
    stringsAsFactors = FALSE
  )
  # block 5 has no entries at all
  err <- capture_sse_error(
    .drawableOmegaBlocks(list(5L), entries, covNames = "om.eta.a")
  )
  expect_s3_class(err, "error")
})

test_that("omegaWishartSpec moment-matches nu from the reported SEs", {
  omega0 <- matrix(c(0.30, 0.05, 0.05, 0.12), 2L, 2L)
  se <- c(0.06, 0.03)
  p <- 2L

  spec <- .omegaWishartSpec(omega0, se)

  expectedNu <- min(p + 3 + 2 * (diag(omega0) / se)^2)
  expect_equal(spec$nu, expectedNu)
  # the spec carries the fitted block; cvPost pre-scaling happens at draw time
  expect_equal(spec$omega0, omega0)
  expect_equal(spec$p, p)
})

test_that("omegaWishartSpec picks the minimum (widest-marginal) nu", {
  omega0 <- diag(c(0.30, 0.12))
  # second element has the relatively larger SE => smaller nu => binding
  spec <- .omegaWishartSpec(omega0, c(0.01, 0.06))
  nuEach <- 2L + 3 + 2 * (diag(omega0) / c(0.01, 0.06))^2
  expect_equal(spec$nu, min(nuEach))
})

test_that("omegaWishartSpec ignores unusable SEs", {
  omega0 <- diag(c(0.30, 0.12))
  # NA SE on the second element => nu comes from the first alone
  spec <- .omegaWishartSpec(omega0, c(0.06, NA_real_))
  expect_equal(spec$nu, 2L + 3 + 2 * (0.30 / 0.06)^2)
})

test_that("omegaWishartSpec returns NULL when no SE is usable", {
  omega0 <- diag(c(0.30, 0.12))
  expect_null(.omegaWishartSpec(omega0, c(NA_real_, NA_real_)))
  expect_null(.omegaWishartSpec(omega0, c(0, 0)))
})

test_that("omegaWishartSpec rejects a degenerate variance", {
  omega0 <- diag(c(0, 0.12))
  err <- capture_sse_error(.omegaWishartSpec(omega0, c(0.06, 0.03)))
  expect_s3_class(err, "error")
})

test_that("omegaWishartSpec names the offending eta in the bad-variance abort", {
  omega0 <- diag(c(0, 0.12))
  dimnames(omega0) <- list(c("eta.ka", "eta.cl"), c("eta.ka", "eta.cl"))
  err <- capture_sse_error(.omegaWishartSpec(omega0, c(0.06, 0.03)))
  expect_s3_class(err, "error")
  expect_match(conditionMessage(err), "eta.ka", fixed = TRUE)
})

test_that("drawOmegaBlock returns a positive-definite symmetric matrix", {
  omega0 <- matrix(c(0.30, 0.05, 0.05, 0.12), 2L, 2L)
  spec <- .omegaWishartSpec(omega0, c(0.06, 0.03))

  set.seed(1234)
  drawn <- .drawOmegaBlock(spec)

  expect_equal(dim(drawn), c(2L, 2L))
  expect_equal(drawn, t(drawn))
  expect_true(all(eigen(drawn, symmetric = TRUE, only.values = TRUE)$values > 0))
})

test_that("drawOmegaBlock recovers the target mean and binding SE", {
  omega0 <- matrix(c(0.30, 0.05, 0.05, 0.12), 2L, 2L)
  targetSe <- c(0.06, 0.03)
  spec <- .omegaWishartSpec(omega0, targetSe)

  set.seed(99)
  draws <- vapply(seq_len(4000L), function(i) diag(.drawOmegaBlock(spec)),
                  numeric(2L))

  # E[Omega] == Omega0 for both elements
  expect_equal(rowMeans(draws), diag(omega0), tolerance = 0.02)

  # the binding (minimum-nu) element recovers its reported SE; here that is
  # element 2, whose SE is relatively largest
  closed <- sqrt(2 * diag(omega0)^2 / (spec$nu - spec$p - 3))
  expect_equal(apply(draws, 1L, stats::sd), closed, tolerance = 0.02)
})

test_that("drawOmega only touches the blocks it is given", {
  omega0 <- diag(c(0.30, 0.12))
  dimnames(omega0) <- list(c("eta.a", "eta.b"), c("eta.a", "eta.b"))

  set.seed(7)
  # Only block 1 is passed, as `.drawableOmegaBlocks()$drawable` would supply.
  # Block 2 is absent from the list and must therefore keep its fitted value.
  # NOTE: do NOT pass an undrawable block here to see what happens -- the
  # contract is that filtering happened upstream, and testing the unfiltered
  # case would encode behaviour callers must never rely on.
  out <- .drawOmega(
    omega0,
    blocks = list(1L),
    se = c(0.06, NA_real_)
  )

  expect_equal(dimnames(out), dimnames(omega0))
  expect_equal(out[2L, 2L], 0.12)
  expect_false(identical(out[1L, 1L], 0.30))
  expect_equal(out[1L, 2L], 0)
})

test_that("drawOmega keeps a drawable block whose SE is unusable", {
  # A block can pass coverage (all entries present and unfixed) yet still have
  # a numerically unusable SE, so no nu can be moment-matched. Coverage and
  # usability are distinct conditions; this block stays at its fitted value.
  omega0 <- matrix(0.30, 1L, 1L, dimnames = list("eta.a", "eta.a"))
  out <- .drawOmega(omega0, blocks = list(1L), se = NA_real_)
  expect_equal(out[1L, 1L], 0.30)
})
