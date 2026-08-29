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
