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
