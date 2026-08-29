# SSE Parameter Uncertainty Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make `parameterSource = "covariance"` draw OMEGA as well as THETA, using the covariance `nlmixr2est` now reports over THETA and OMEGA jointly — fixing a currently broken code path.

**Architecture:** Two selectable draw modes sharing one coverage policy.
`covarianceDraw = "independent_iw"` (**the default**) draws THETA
multivariate-Normal from the theta sub-block and OMEGA from a mean-centred
inverse-Wishart per block, independently, delegating the OMEGA draw to
`rxode2::cvPost()`. `covarianceDraw = "joint"` (**opt-in**) draws THETA and
OMEGA together from the full covariance on a log-Cholesky-transformed scale,
incorporating the THETA↔OMEGA covariance through a first-order delta
approximation. Both are positive-definite by construction.

**Neither mode is NONMEM `$PRIOR NWPRI`, and neither claims PsN parity.** An
earlier revision of this plan called the independent mode `"nwpri"` and
asserted parity; that was wrong (rxode2 documents that NWPRI does not reduce
to a textbook inverse-Wishart), and it has been renamed.

**Tech Stack:** R package, `rxode2::cvPost` for the inverse-Wishart draw, `nlmixr2est` covariance output, testthat 3e, checkmate, cli.

---

## Background For The Implementer

You need no prior context. Read this section fully — it contains verified facts you must not re-derive.

### The defect being fixed

`nlmixr2est` now CAN report a **joint** fixed/random-effect covariance in `fit$cov` -- for FOCEI `"r"`, `"s"`, `"r,s"`, `"analytic"` and SAEM `"fim"`, `"sa"` -- whenever the full-covariance calculation and installation succeed. It is not guaranteed: see the coverage policy below. OMEGA entries appear alongside the thetas:

```
dimnames(fit$cov)[[1]]
#> "tka" "tcl" "tv" "add.sd" "om.eta.ka" "cov.eta.cl.eta.ka" "om.eta.cl" "om.eta.v"
```

`.alignedCovariance()` assumes `fit$cov` is theta-only and aborts on the first name absent from `fit$theta`:

```
Error in .alignedCovariance(fit) :
  `fit$cov` contains theta name "om.eta.ka" that are not present in `fit$theta`.
```

This reproduces with the exact model and control settings in the package README, so covariance-mode SSE currently fails for ordinary usage.

### Verified facts — do not re-derive these

**1. OMEGA naming in `fit$cov`.** Produced by an unexported `nlmixr2est` helper. The rule is:

- diagonal (`row == col`) → `paste0("om.", rowName)`
- off-diagonal → `paste0("cov.", rowName, ".", colName)`

Confirmed against a live fit with a 2×2 block plus a separate diagonal eta.

**2. `.uiOmegaInfo()` gives you the mapping.** It already exists in `R/sse-helpers.R`. **It must be passed `fit$ui`, NOT `fit`** — passing a fit returns zero rows, silently. Verified output:

```
  row col rowName colName   fix
1   1   1  eta.ka  eta.ka FALSE
2   2   1  eta.cl  eta.ka FALSE
3   2   2  eta.cl  eta.cl FALSE
4   3   3   eta.v   eta.v FALSE
```

Reuse it. Do not write a new iniDf parser.

**3. `rxode2::cvPost()` is exported** and implements the inverse-Wishart draw. Do not use `stats::rWishart` directly.

**4. The inverse-Wishart moment-match formula.** For `Omega ~ InvWishart(Psi, nu)` over a `p × p` block, setting `Psi = (nu - p - 1) * Omega0` gives `E[Omega] = Omega0` exactly and:

```
Var(Omega_ii) = 2 * Omega0_ii^2 / (nu - p - 3)
```

Solving against the reported standard error `SE_i`:

```
nu_i = p + 3 + 2 * (Omega0_ii / SE_i)^2
```

Take `nu = min(nu_i)` across the block's diagonal elements (the widest choice for the constructed diagonal marginals; inverse-Wishart has one `nu` per matrix).

**Verified numerically.** 2×2 block, `Omega0 = [[0.30, 0.05], [0.05, 0.12]]`, target SEs `(0.06, 0.03)`: 200,000 draws gave `E[Omega_ii] = (0.300, 0.120)` against target `(0.300, 0.120)`, and `SD = (0.0751, 0.0300)` against closed-form `(0.0750, 0.0300)`.

**Known, accepted limitation:** only the binding (minimum-`nu`) element matches its SE exactly. In that verification the non-binding element got `SD = 0.075` against its `0.06` target — 25% wider. `min(nu_i)` always errs wider, never narrower, for the CONSTRUCTED DIAGONAL MARGINALS. That is not a general conservatism guarantee -- it says nothing about non-linear predictions, correlations, delta-OFV, Type-I error, or power. This is by design; document it, do not try to fix it.

**5. Drawing — use `rxode2::cvPost()`, do not hand-roll it.** It is exported
and already implements the inverse-Wishart (`type = "invWishart"`, plus
`"lkj"` and `"separation"`).

Its convention is `Psi = nu * omega`, so it is **not** mean-centred:
`E[Omega*] = nu/(nu - p - 1) * omega`. Verified empirically — at `nu = 20,
p = 2` the observed ratio is 1.179 against a predicted 1.176; at `nu = 200`,
1.014 against 1.015. Pre-scale to cancel it:

```r
scaled <- omega * (nu - p - 1) / nu
omegaBlock <- rxode2::cvPost(nu, scaled, n = 1L)
```

Verified: this recovers `E[Omega] = 0.2998/0.1200` against a target of
`0.30/0.12`, an SD matching the closed form `sqrt(2*omega_ii^2/(nu-p-3))`, and
positive-definite draws throughout.

**6. Residual error needs no special handling.** In nlmixr2 residual-error parameters (`add.sd`, `prop.sd`, …) are ordinary **thetas** and appear in the theta block of `fit$cov` with full SEs. They are drawn by the THETA step automatically. `fit$sigma` returns `named numeric(0)` — that is a legacy accessor artifact, NOT missing uncertainty. Do not add a SIGMA draw.

### Design decisions already settled — do NOT redesign

- **Two modes, `"independent_iw"` the default.** `"joint"` is opt-in because
  its non-linear back-transform inflates OMEGA means by `exp(SE²/2Ω̂²)` — 2% at
  20% relative SE but **64.9% at 100%**, worst exactly where uncertainty
  propagation matters most. Do not make `"joint"` the default.
- **Delegate the inverse-Wishart to `rxode2::cvPost()`.** Do not hand-roll it
  with `rWishart`. Pre-scale as described in Task 4 to cancel `cvPost`'s
  `Psi = nu * omega` convention.
- **Coverage policy** (below) is shared by both modes and is not negotiable:
  partial coverage is a supported degraded mode, and any block with a fixed or
  uncovered element is held entirely fixed.
- **Per-block `nu`**, not one global `nu`.
- **No backward compatibility of draw values.** Same-seed results will differ
  from previous releases; that is expected and gets a NEWS entry.

### Coverage policy — read before writing any draw code

`fit$cov` does not always cover every parameter. Decide coverage **per
parameter and per block**, never all-or-nothing:

- **Always draw covered thetas**, even when no OMEGA block is drawable. A
  theta-only `fit$cov` is a *supported configuration*, not a failure —
  `foceiControl(covFull = FALSE)` requests it, and the installer falls back to
  it when full-covariance components are missing or non-positive-definite.
- **Draw only fully covered OMEGA blocks.**
- **Hold an entire block fixed if any of its declared elements is fixed or
  missing from `fit$cov`.** `nlmixr2est` drops fixed OMEGA elements from
  `fit$cov` (`.foceiOmegaPairs()` filters on `.omegaFixed()`) while the
  declared topology still lists them, so a correlated block with one fixed
  component *will* have missing entries. Do **not** draw the free sub-matrix
  and splice fitted values back in — replacing entries of a positive-definite
  draw can yield a non-positive-definite matrix.
- **Record everything not drawn** in `parameterPartition$fixed`.

### A silent R trap you must avoid

`x[-seq_len(n)]` returns an **empty** vector when `n == 0`, because
`-integer(0)` is `integer(0)`. Verified: `c(10,20,30)[-seq_len(0)]` has length
0, not 3. Any code splitting a joint draw into theta and omega parts by
negative indexing will silently discard the entire OMEGA draw for an
OMEGA-only covariance. Split by explicit positive indices instead, and test the
zero-theta case.

### Current code

- `R/sse-helpers.R:363` `.validateCovarianceFit()` — checks `fit$cov` is a non-empty square PD matrix. Unchanged by this plan.
- `R/sse-helpers.R:399` `.paramSetFromFit()` — returns `list(theta, omega, sigma)`.
- `R/sse-helpers.R:799` `.alignedCovariance()` — the function that aborts. Task 5 rewrites it.
- `R/sse-helpers.R:839` `.resolveCovarianceParameterSets()` — does the per-replicate draw. Task 6 rewrites it.
- `R/sse-helpers.R:1115` `.uiOmegaInfo()` — reuse as-is.

Line numbers are from the current branch tip and may drift as you work; always locate by function name.

## File Structure

- `R/sse-omega-draw.R` — **new file.** Shared OMEGA plumbing (name construction, block detection, raw-results labels) plus the `"independent_iw"` inverse-Wishart machinery.
- `R/sse-omega-joint.R` — **new file.** The `"joint"` mode: log-Cholesky transform, numeric Jacobian, joint draw spec and draw. Kept separate from the Wishart code because the two modes are independent implementations of the same contract, and each is easier to reason about alone.
- `R/sse-helpers.R` — `.alignedCovariance()` rewritten to partition rather than abort; `.resolveCovarianceParameterSets()` rewritten to dispatch on `covarianceDraw`.
- `R/sse-control.R` — new `covarianceDraw` argument.
- `tests/testthat/test-sse-omega-draw.R` — **new file.** Shared plumbing + `"independent_iw"` tests.
- `tests/testthat/test-sse-omega-joint.R` — **new file.** `"joint"` mode tests.
- `tests/testthat/test-control.R` — `covarianceDraw` validation.
- `NEWS.md`, `README.md`, `man/runSSEControl.Rd` — documentation.

`R/sse-helpers.R` is already ~2000 lines, so both new concerns go in their own
files rather than growing it further.

## Conventions

- Tests begin with `skip_on_cran()`. `tests/testthat/test-run-sse.R` and `test-control.R` already have it; a NEW test file needs one.
- Errors use the package-local `.abortSSE()`, which takes a cli-formatted string.
- `capture_sse_error(expr)` (in `tests/testthat/helper-fake-fit.R`) returns the condition or `NULL`.
- Run one file: `NOT_CRAN=true Rscript -e "devtools::load_all(quiet=TRUE); testthat::test_file('tests/testthat/<file>', reporter='summary')"`
- **Do NOT run `devtools::document()`** — this environment has roxygen2 8.1.0 against a package pinned to `RoxygenNote: 7.3.3`, and it rewrites many unrelated files. Use `@noRd` on all new internal functions so no `.Rd` is needed.

## Pre-Existing Failures

`devtools::test()` fails before this plan starts. Use this list to tell your regressions from inherited breakage. **Do not fix these.**

1. `test-run-sse.R` — `simulationRecord preserves row identifiers…` → `bad matrix specification`
2. `test-run-sse.R` — `runSSE resumes cached partial runs…` → `numbers of columns of arguments do not match`
3. `test-run-sse.R` — `runSSE addModels reuses saved datasets…` → `Could not recover the original estimation data…`
4. `test-run-sse.R` — `runSSE supports raw-results and covariance parameter sources end to end` → **SKIPS**, masked by a `cli` `.envir` defect in `nlmixr2utils:::.abortRawResults()`
5. `test-sse-power-vignette-data.R` — two failures asserting README mentions `vignette("runSSE"…)` / `sse-power`

Baseline at plan start: `[ FAIL 5 | WARN 0 | SKIP 1 | PASS 137 ]`.

Note failure #4 is the integration test for the very mode this plan fixes. It is masked by a *sibling package* bug and will likely still skip. Do not chase it.

### Assertions and documentation that this work makes false

Do **not** assume the expected-failure list is unchanged. Once OMEGA is drawn:

| Location | What happens |
| --- | --- |
| `tests/testthat/test-run-sse.R:888` — asserts drawn OMEGA has exactly one unique value | **Now false**, but *currently masked*: it sits inside the end-to-end test that `skip()`s (failure #4). Task 12 corrects it anyway — left alone it becomes a trap the moment that sibling-package defect is fixed. |
| `tests/testthat/test-run-sse.R:168` — "draws theta and keeps omega fixed" | **Stays correct.** Its `fake_sse_fit()` fixture has a theta-only `cov` (`tka`, `tcl` only), so the coverage policy holds OMEGA fixed. Keep it — it is now the theta-only regression test. Do not "fix" it. |
| `R/run-sse.R:291` — cli message "OMEGA and SIGMA stay at the fitted point estimates" | Now false. Task 12. |
| `README.md:97` — "OMEGA and SIGMA stay fixed" | Now false. Task 12. |
| `vignettes/runSSE.Rmd:203, 219` — "theta-only uncertainty" / "For full uncertainty in OMEGA and SIGMA as well, use …" | Now false. Task 12. |

Line numbers are from the current branch tip and may drift; locate by content.

---

### Task 1: OMEGA covariance-name construction

**Files:**
- Create: `R/sse-omega-draw.R`
- Create: `tests/testthat/test-sse-omega-draw.R`

- [ ] **Step 1: Write the failing test**

Create `tests/testthat/test-sse-omega-draw.R`:

```r
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
```

- [ ] **Step 2: Run to verify it fails**

Run: `NOT_CRAN=true Rscript -e "devtools::load_all(quiet=TRUE); testthat::test_file('tests/testthat/test-sse-omega-draw.R', reporter='summary')"`
Expected: FAIL — `could not find function ".omegaCovName"`.

- [ ] **Step 3: Implement**

Create `R/sse-omega-draw.R`:

```r
# Inverse-Wishart machinery for drawing OMEGA with uncertainty.
#
# Shares NWPRI's broad independence factorization -- THETA and OMEGA drawn from INDEPENDENT
# distributions (Normal and inverse-Wishart respectively), so the THETA<->OMEGA
# cross-terms present in fit$cov are deliberately not used here.

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
```

- [ ] **Step 4: Run to verify it passes**

Same command as Step 2. Expected: PASS.

- [ ] **Step 5: Verify against a real fit**

This confirms the naming rule really matches `nlmixr2est`, which is the assumption the whole feature rests on.

```bash
NOT_CRAN=true Rscript -e '
devtools::load_all(quiet = TRUE)
suppressMessages(library(nlmixr2))
m <- function() {
  ini({
    tka <- log(1.57); tcl <- log(2.72); tv <- log(31.5)
    eta.ka + eta.cl ~ c(0.6, 0.1, 0.3)
    eta.v ~ 0.2
    add.sd <- 0.7
  })
  model({
    ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
    cp <- linCmt(); cp ~ add(add.sd)
  })
}
fit <- suppressMessages(suppressWarnings(nlmixr2(m, nlmixr2data::theo_sd,
  est = "focei", control = list(print = 0L, covMethod = "r", eval.max = 40L))))
tab <- nlmixr2sse:::.omegaEntryTable(nlmixr2sse:::.uiOmegaInfo(fit$ui))
covNames <- dimnames(fit$cov)[[1]]
cat("constructed:", paste(tab$covName, collapse = ", "), "\n")
cat("in fit$cov :", paste(intersect(tab$covName, covNames), collapse = ", "), "\n")
stopifnot(all(tab$covName %in% covNames))
cat("ALL CONSTRUCTED NAMES PRESENT IN fit$cov\n")
'
```

Expected: ends with `ALL CONSTRUCTED NAMES PRESENT IN fit$cov`.

If this fails, STOP and report — the naming convention has changed in `nlmixr2est` and the whole approach needs revisiting.

- [ ] **Step 6: Commit**

```bash
git add R/sse-omega-draw.R tests/testthat/test-sse-omega-draw.R
git commit -m "feat: construct nlmixr2est-style OMEGA covariance names"
```

---

### Task 2: OMEGA block detection and the coverage policy

**Files:**
- Modify: `R/sse-omega-draw.R`
- Modify: `tests/testthat/test-sse-omega-draw.R`

An inverse-Wishart draw applies to one OMEGA block. Blocks are the connected components of the off-diagonal structure: two etas share a block if an off-diagonal entry links them (directly or transitively).

- [ ] **Step 1: Write the failing test**

Append to `tests/testthat/test-sse-omega-draw.R`:

```r
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
```

- [ ] **Step 2: Run to verify it fails**

Same command as Task 1 Step 2. Expected: FAIL — `could not find function ".omegaBlocks"`.

- [ ] **Step 3: Implement**

Append to `R/sse-omega-draw.R`:

```r
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
  split(seq_len(nEta), factor(roots, levels = unique(roots)))
}
```

Note `split()` returns a named list; the tests compare against an unnamed list, so strip names:

```r
  unname(split(seq_len(nEta), factor(roots, levels = unique(roots))))
```

Use that as the final line instead.

- [ ] **Step 4: Run to verify it passes**

Same command. Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add R/sse-omega-draw.R tests/testthat/test-sse-omega-draw.R
git commit -m "feat: detect independent OMEGA blocks for Wishart draws"
```

- [ ] **Step 6: Write the failing coverage-policy test**

This implements the coverage policy from the preamble. **Both modes consume
this one helper** — neither may re-derive drawability from standard errors.

Append to `tests/testthat/test-sse-omega-draw.R`:

```r
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
```

- [ ] **Step 7: Run to verify it fails**

Expected: FAIL — `could not find function ".drawableOmegaBlocks"`.

- [ ] **Step 8: Implement the coverage policy**

Append to `R/sse-omega-draw.R`:

```r
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
```

- [ ] **Step 9: Run to verify it passes**

Expected: PASS.

- [ ] **Step 10: Commit**

```bash
git add R/sse-omega-draw.R tests/testthat/test-sse-omega-draw.R
git commit -m "feat: hold OMEGA blocks fixed unless fully covered by fit\$cov"
```

---

### Task 3 (`"independent_iw"` mode): Moment-matched inverse-Wishart spec per block

**Files:**
- Modify: `R/sse-omega-draw.R`
- Modify: `tests/testthat/test-sse-omega-draw.R`

- [ ] **Step 1: Write the failing test**

Append to `tests/testthat/test-sse-omega-draw.R`:

```r
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
```

- [ ] **Step 2: Run to verify it fails**

Expected: FAIL — `could not find function ".omegaWishartSpec"`.

- [ ] **Step 3: Implement**

Append to `R/sse-omega-draw.R`:

```r
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
    .abortSSE(
      paste0(
        "OMEGA variance {.val {bad}} is zero or non-finite, so its uncertainty ",
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
```

- [ ] **Step 4: Run to verify it passes**

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add R/sse-omega-draw.R tests/testthat/test-sse-omega-draw.R
git commit -m "feat: moment-match inverse-Wishart df from OMEGA standard errors"
```

---

### Task 4 (`"independent_iw"` mode): Draw OMEGA

**Files:**
- Modify: `R/sse-omega-draw.R`
- Modify: `tests/testthat/test-sse-omega-draw.R`

- [ ] **Step 1: Write the failing test**

Append to `tests/testthat/test-sse-omega-draw.R`:

```r
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
```

- [ ] **Step 2: Run to verify it fails**

Expected: FAIL — `could not find function ".drawOmegaBlock"`.

- [ ] **Step 3: Implement**

Append to `R/sse-omega-draw.R`:

```r
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
```

- [ ] **Step 4: Run to verify it passes**

Expected: PASS. The statistical tests use 4000 draws with a 0.02 tolerance; if one is marginally flaky, re-run once to confirm before changing anything — do NOT loosen the tolerance to hide a real problem.

- [ ] **Step 5: Commit**

```bash
git add R/sse-omega-draw.R tests/testthat/test-sse-omega-draw.R
git commit -m "feat: draw OMEGA blocks from moment-matched inverse-Wisharts"
```

---

### Task 5: Partition `fit$cov` instead of aborting

**Files:**
- Modify: `R/sse-helpers.R` (`.alignedCovariance()`)
- Modify: `tests/testthat/test-sse-omega-draw.R`

This is the task that fixes the reported defect.

- [ ] **Step 1: Write the failing test**

Append to `tests/testthat/test-sse-omega-draw.R`:

```r
test_that("alignedCovariance partitions theta and omega entries", {
  nm <- c("tka", "tcl", "om.eta.ka", "cov.eta.cl.eta.ka", "om.eta.cl")
  cov <- diag(c(0.04, 0.09, 0.01, 0.002, 0.005))
  dimnames(cov) <- list(nm, nm)

  omega <- matrix(c(0.3, 0.02, 0.02, 0.2), 2L, 2L)
  dimnames(omega) <- list(c("eta.ka", "eta.cl"), c("eta.ka", "eta.cl"))

  fit <- list(
    theta = c(tka = 0.45, tcl = 1.0),
    omega = omega,
    cov = cov,
    ui = list(iniDf = data.frame(
      name = c("tka", "tcl", "eta.ka", "(eta.ka,eta.cl)", "eta.cl"),
      ntheta = c(1L, 2L, NA, NA, NA),
      neta1 = c(NA, NA, 1L, 2L, 2L),
      neta2 = c(NA, NA, 1L, 1L, 2L),
      fix = rep(FALSE, 5L),
      stringsAsFactors = FALSE
    ))
  )

  aligned <- .alignedCovariance(fit)

  expect_equal(aligned$drawNames, c("tka", "tcl"))
  expect_equal(dim(aligned$cov), c(2L, 2L))
  expect_equal(aligned$omegaSe[["eta.ka"]], sqrt(0.01))
  expect_equal(aligned$omegaSe[["eta.cl"]], sqrt(0.005))
})

test_that("alignedCovariance no longer aborts on omega names", {
  nm <- c("tka", "om.eta.ka")
  cov <- diag(c(0.04, 0.01))
  dimnames(cov) <- list(nm, nm)
  omega <- matrix(0.3, 1L, 1L, dimnames = list("eta.ka", "eta.ka"))

  fit <- list(
    theta = c(tka = 0.45),
    omega = omega,
    cov = cov,
    ui = list(iniDf = data.frame(
      name = c("tka", "eta.ka"),
      ntheta = c(1L, NA),
      neta1 = c(NA, 1L),
      neta2 = c(NA, 1L),
      fix = c(FALSE, FALSE),
      stringsAsFactors = FALSE
    ))
  )

  expect_silent(aligned <- .alignedCovariance(fit))
  expect_equal(aligned$drawNames, "tka")
})

test_that("alignedCovariance still rejects a genuinely unknown name", {
  nm <- c("tka", "mystery")
  cov <- diag(c(0.04, 0.01))
  dimnames(cov) <- list(nm, nm)

  fit <- list(
    theta = c(tka = 0.45),
    omega = matrix(numeric(0), 0L, 0L),
    cov = cov,
    ui = list(iniDf = data.frame(
      name = "tka", ntheta = 1L, neta1 = NA_integer_, neta2 = NA_integer_,
      fix = FALSE, stringsAsFactors = FALSE
    ))
  )

  err <- capture_sse_error(.alignedCovariance(fit))
  expect_s3_class(err, "error")
  expect_match(conditionMessage(err), "mystery")
})
```

- [ ] **Step 2: Run to verify it fails**

Expected: FAIL — the current `.alignedCovariance()` aborts on `om.eta.ka`, and returns no `omegaSe`.

- [ ] **Step 3: Implement**

In `R/sse-helpers.R`, replace the body of `.alignedCovariance()` from the `missing_theta` check to the end of the function. The current tail is:

```r
  missing_theta <- setdiff(cov_names, names(theta))
  if (length(missing_theta) > 0L) {
    .abortSSE(
      "{.arg fit$cov} contains theta name{?s} {.val {missing_theta}} that are not present in {.arg fit$theta}."
    )
  }

  list(
    theta = theta,
    cov = cov_mat,
    drawNames = cov_names
  )
}
```

Replace with:

```r
  # fit$cov now carries OMEGA entries (om.<eta> / cov.<eta>.<eta>) alongside the
  # thetas. Partition them: thetas are drawn multivariate-Normal from their own
  # sub-block, OMEGA is drawn separately from an inverse-Wishart, and the
  # THETA<->OMEGA cross-terms are deliberately discarded (this mode shares NWPRI's independence factorization, though it is not
  # NWPRI -- the OMEGA density differs).
  omega_entries <- .omegaEntryTable(.uiOmegaInfo(fit[["ui"]]))
  omega_names <- omega_entries$covName

  theta_names <- intersect(cov_names, names(theta))
  unknown <- setdiff(cov_names, c(theta_names, omega_names))
  if (length(unknown) > 0L) {
    .abortSSE(
      "{.arg fit$cov} contains name{?s} {.val {unknown}} that match neither {.arg fit$theta} nor an OMEGA entry."
    )
  }

  # per-eta standard errors, from the diagonal OMEGA entries present in fit$cov
  omega_mat <- .fitField(fit, "omega")
  eta_names <- if (is.matrix(omega_mat)) rownames(omega_mat) else character(0)
  omega_se <- stats::setNames(
    rep(NA_real_, length(eta_names)),
    eta_names
  )
  diag_entries <- omega_entries[omega_entries$diagonal, , drop = FALSE]
  for (i in seq_len(nrow(diag_entries))) {
    nm <- diag_entries$covName[[i]]
    eta <- diag_entries$rowName[[i]]
    if (isTRUE(diag_entries$fix[[i]])) {
      next
    }
    if (nm %in% cov_names && eta %in% eta_names) {
      omega_se[[eta]] <- sqrt(cov_mat[nm, nm])
    }
  }

  list(
    theta = theta,
    cov = cov_mat[theta_names, theta_names, drop = FALSE],
    # the whole matrix, kept so the "joint" mode can pull a theta+omega
    # sub-block including the cross-terms that "independent_iw" discards
    fullCov = cov_mat,
    drawNames = theta_names,
    omegaSe = omega_se,
    omegaEntries = omega_entries
  )
}
```

- [ ] **Step 4: Run to verify it passes**

Expected: PASS.

- [ ] **Step 5: Run the full suite**

Run: `NOT_CRAN=true Rscript -e "devtools::load_all(quiet=TRUE); devtools::test()"`
Expected: FAIL still 5 (the pre-existing list), PASS higher. Anything else is your regression.

- [ ] **Step 6: Commit**

```bash
git add R/sse-helpers.R tests/testthat/test-sse-omega-draw.R
git commit -m "fix: partition theta and omega entries in fit\$cov instead of aborting"
```

---

### Task 6 (`"joint"` mode): Log-Cholesky transform

**Files:**
- Create: `R/sse-omega-joint.R`
- Create: `tests/testthat/test-sse-omega-joint.R`

The opt-in `"joint"` mode draws THETA and OMEGA together. To keep every draw
positive-definite, OMEGA is mapped to an unconstrained log-Cholesky vector
first: `L = t(chol(Omega))`, take `L`'s lower-triangular elements in
column-major order, and replace the diagonal entries with their logs. Inverting
exponentiates the diagonal, so `L %*% t(L)` is PD for **any** input vector.

- [ ] **Step 1: Write the failing test**

Create `tests/testthat/test-sse-omega-joint.R`:

```r
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
```

- [ ] **Step 2: Run to verify it fails**

Run: `NOT_CRAN=true Rscript -e "devtools::load_all(quiet=TRUE); testthat::test_file('tests/testthat/test-sse-omega-joint.R', reporter='summary')"`
Expected: FAIL — `could not find function ".omegaToPhi"`.

- [ ] **Step 3: Implement**

Create `R/sse-omega-joint.R`:

```r
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
```

- [ ] **Step 4: Run to verify it passes**

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add R/sse-omega-joint.R tests/testthat/test-sse-omega-joint.R
git commit -m "feat: log-Cholesky transform for positive-definite OMEGA draws"
```

---

### Task 7 (`"joint"` mode): Delta-method Jacobian

**Files:**
- Modify: `R/sse-omega-joint.R`
- Modify: `tests/testthat/test-sse-omega-joint.R`

`fit$cov` describes uncertainty on the **natural** OMEGA scale. To draw on the
transformed scale it must be carried across with the delta method, which needs
`J = d(phi)/d(omega)`. The number of free OMEGA elements is small, so central
finite differences are accurate and cheap, and avoid a new dependency.

- [ ] **Step 1: Write the failing test**

Append to `tests/testthat/test-sse-omega-joint.R`:

```r
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
```

- [ ] **Step 2: Run to verify it fails**

Expected: FAIL — `could not find function ".numericJacobian"`.

- [ ] **Step 3: Implement**

Append to `R/sse-omega-joint.R`:

```r
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
```

Callers must supply `context` so the error names the block rather than an
opaque element number. In `.jointDrawSpec()`, change the Jacobian call to:

```r
    jb <- .numericJacobian(
      function(w) .omegaToPhi(.vecToOmega(w, p)),
      .omegaToVec(b$omega),
      context = paste(rownames(b$omega), collapse = ", ")
    )
```

- [ ] **Step 4: Run to verify it passes**

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add R/sse-omega-joint.R tests/testthat/test-sse-omega-joint.R
git commit -m "feat: central-difference Jacobian for delta-method transform"
```

---

### Task 8 (`"joint"` mode): Build the joint spec and draw

**Files:**
- Modify: `R/sse-omega-joint.R`
- Modify: `tests/testthat/test-sse-omega-joint.R`

Assemble `Sigma_T = B Sigma B'` with `B = blockdiag(I, J)`, factor it once, and
draw `(theta, phi)` together per replicate.

- [ ] **Step 1: Write the failing test**

Append to `tests/testthat/test-sse-omega-joint.R`:

```r
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
  expect_equal(stats::sd(om11), 0.060, tolerance = 0.01)
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
```

- [ ] **Step 2: Run to verify it fails**

Expected: FAIL — `could not find function ".jointDrawSpec"`.

- [ ] **Step 3: Implement**

Append to `R/sse-omega-joint.R`:

```r
#' Build the transformed joint draw specification
#'
#' `blocks` is a list of `list(omega = <block matrix>, index = <eta indices>)`.
#' `sigma` is the covariance over `c(theta, <block omega vectors in order>)` on
#' the NATURAL scale, already subset to the drawn entries.
#'
#' @return list(thetaNames, thetaMean, phiMean, blocks, chol, nTheta)
#' @noRd
.jointDrawSpec <- function(theta, blocks, sigma) {
  nTheta <- length(theta)
  phiMean <- unlist(lapply(blocks, function(b) .omegaToPhi(b$omega)),
                    use.names = FALSE)

  # Jacobian of the stacked natural->transformed map, block-diagonal because
  # blocks are structurally independent of one another.
  jac <- NULL
  for (b in blocks) {
    p <- nrow(b$omega)
    jb <- .numericJacobian(
      function(w) .omegaToPhi(.vecToOmega(w, p)),
      .omegaToVec(b$omega)
    )
    jac <- if (is.null(jac)) {
      jb
    } else {
      rbind(
        cbind(jac, matrix(0, nrow(jac), ncol(jb))),
        cbind(matrix(0, nrow(jb), ncol(jac)), jb)
      )
    }
  }

  bMat <- if (is.null(jac)) {
    diag(nTheta)
  } else {
    rbind(
      cbind(diag(nTheta), matrix(0, nTheta, ncol(jac))),
      cbind(matrix(0, nrow(jac), nTheta), jac)
    )
  }

  sigmaT <- bMat %*% sigma %*% t(bMat)
  sigmaT <- 0.5 * (sigmaT + t(sigmaT))

  ev <- suppressWarnings(
    eigen(sigmaT, symmetric = TRUE, only.values = TRUE)$values
  )
  if (!all(is.finite(ev)) || min(ev) <= 0) {
    .abortSSE(
      paste0(
        "The delta-method transformed covariance is not positive-definite, so ",
        "a joint THETA/OMEGA draw is not possible. Use ",
        "{.code runSSEControl(covarianceDraw = \"independent_iw\")}, or ",
        "{.code parameterSource = \"fixed\"}."
      )
    )
  }

  list(
    thetaNames = names(theta),
    thetaMean = unname(theta),
    phiMean = phiMean,
    blocks = blocks,
    chol = chol(sigmaT),
    nTheta = nTheta
  )
}

#' Take one joint draw and back-transform it
#'
#' @return list(theta = named numeric, omega = list of block matrices)
#' @noRd
.drawJoint <- function(spec) {
  mu <- c(spec$thetaMean, spec$phiMean)
  drawn <- mu + as.numeric(stats::rnorm(length(mu)) %*% spec$chol)

  # NOTE: split by explicit positive indices. `drawn[-seq_len(nTheta)]` returns
  # an EMPTY vector when nTheta == 0 (because -integer(0) is integer(0)), which
  # would silently discard the entire OMEGA draw for an OMEGA-only covariance.
  theta <- if (spec$nTheta > 0L) {
    stats::setNames(drawn[seq_len(spec$nTheta)], spec$thetaNames)
  } else {
    stats::setNames(numeric(0), character(0))
  }

  phi <- if (length(drawn) > spec$nTheta) {
    drawn[(spec$nTheta + 1L):length(drawn)]
  } else {
    numeric(0)
  }
  omega <- list()
  at <- 0L
  for (b in spec$blocks) {
    p <- nrow(b$omega)
    n <- p * (p + 1L) / 2L
    omega[[length(omega) + 1L]] <- .phiToOmega(phi[at + seq_len(n)], p)
    at <- at + n
  }

  list(theta = theta, omega = omega)
}
```

- [ ] **Step 4: Run to verify it passes**

Expected: PASS. These are Monte-Carlo tests with generous tolerances; if one is
marginally flaky, re-run once before changing anything — do NOT loosen a
tolerance to hide a real problem.

- [ ] **Step 5: Commit**

```bash
git add R/sse-omega-joint.R tests/testthat/test-sse-omega-joint.R
git commit -m "feat: joint THETA/OMEGA draw incorporating estimated covariance"
```

---

### Task 9: Add the `covarianceDraw` control argument

**Files:**
- Modify: `R/sse-control.R`
- Modify: `tests/testthat/test-control.R`

- [ ] **Step 1: Write the failing test**

Append to `tests/testthat/test-control.R`:

```r
test_that("runSSEControl exposes covarianceDraw defaulting to independent_iw", {
  expect_equal(runSSEControl()$covarianceDraw, "independent_iw")
  expect_equal(
    runSSEControl(
      parameterSource = "covariance",
      covarianceDraw = "joint"
    )$covarianceDraw,
    "joint"
  )
})

test_that("runSSEControl exposes omegaRseWarn", {
  expect_equal(runSSEControl()$omegaRseWarn, 0.5)
  expect_equal(
    runSSEControl(
      parameterSource = "covariance",
      omegaRseWarn = 0.25
    )$omegaRseWarn,
    0.25
  )
  err <- capture_sse_error(
    runSSEControl(parameterSource = "covariance", omegaRseWarn = -1)
  )
  expect_s3_class(err, "error")
})

test_that("runSSEControl rejects an unknown covarianceDraw", {
  err <- capture_sse_error(
    runSSEControl(parameterSource = "covariance", covarianceDraw = "wishart")
  )
  expect_s3_class(err, "error")
})

test_that("covarianceDraw requires parameterSource = covariance", {
  err <- capture_sse_error(
    runSSEControl(parameterSource = "fixed", covarianceDraw = "independent_iw")
  )
  expect_s3_class(err, "error")
  expect_match(conditionMessage(err), "covarianceDraw")
})
```

- [ ] **Step 2: Run to verify it fails**

Run: `NOT_CRAN=true Rscript -e "devtools::load_all(quiet=TRUE); testthat::test_file('tests/testthat/test-control.R', reporter='summary')"`
Expected: FAIL — `covarianceDraw` is `NULL`.

- [ ] **Step 3: Implement**

In `R/sse-control.R`, add the formal immediately after `parameterSource`:

```r
  parameterSource = c("fixed", "rawres", "covariance"),
  covarianceDraw = c("independent_iw", "joint"),
  omegaRseWarn = 0.5,
```

Note the order: `match.arg()` takes the **first** element as the default, so
`"independent_iw"` must come first. `"joint"` is opt-in because its
back-transform inflates OMEGA means by `exp(SE²/2Ω̂²)` — 64.9% at 100%
relative SE.
```

Immediately after the existing `parameterSource <- match.arg(parameterSource)`, add:

```r
  covarianceDrawMissing <- missing(covarianceDraw)
  covarianceDraw <- match.arg(covarianceDraw)
```

Add this gate alongside the other `parameterSource` gates (next to the
`inFilter` one):

```r
  if (!covarianceDrawMissing && parameterSource != "covariance") {
    .abortSSE(
      "{.arg covarianceDraw} requires {.arg parameterSource = \"covariance\"}."
    )
  }
```

Add to the returned `structure()` list, after `parameterSource = parameterSource,`:

```r
      covarianceDraw = covarianceDraw,
```

Add the roxygen `@param` after the `@param parameterSource` block:

```r
#' @param covarianceDraw How `parameterSource = "covariance"` draws parameters.
#'   `"independent_iw"` (the default) draws THETA multivariate-Normal and OMEGA
#'   from a mean-centred inverse-Wishart per OMEGA block, independently of each
#'   other. `"joint"` instead draws THETA and OMEGA together from `fit$cov` on a
#'   log-Cholesky-transformed scale, incorporating the estimated THETA/OMEGA
#'   covariance through a first-order delta approximation. Both give
#'   positive-definite OMEGA draws.
#'
#'   `"joint"` uses information `"independent_iw"` discards, but its non-linear
#'   back-transform inflates OMEGA means by roughly `exp(SE^2 / (2 * Omega^2))`
#'   — about 2% at 20% relative standard error, but 65% at 100% — so it is
#'   opt-in rather than the default. Neither mode implements NONMEM's
#'   `$PRIOR NWPRI`, and neither claims PsN parity.
#' @param omegaRseWarn Relative-standard-error threshold above which a drawn
#'   OMEGA variance triggers a weak-identification warning naming the block and
#'   eta. Defaults to `0.5`.
```

- [ ] **Step 4: Run to verify it passes**

Expected: PASS.

- [ ] **Step 5: Hand-patch the documentation**

Do NOT run `devtools::document()`. Add the entries to `man/runSSEControl.Rd`
by hand. **Both** new arguments must be patched, not just `covarianceDraw`:

To the `\usage{}` block, in argument order:

```
  covarianceDraw = c("independent_iw", "joint"),
  omegaRseWarn = 0.5,
```

To `\arguments{}`, an `\item{covarianceDraw}{...}` **and** an
`\item{omegaRseWarn}{...}`, each matching its roxygen text above.

A `\usage{}` block that omits an argument present in the function signature is
an `R CMD check` WARNING ("Undocumented arguments"), so missing `omegaRseWarn`
here would surface later as a check failure rather than a silent gap.

Verify it still parses:

```bash
Rscript -e 'invisible(tools::parse_Rd("man/runSSEControl.Rd")); cat("Rd parses OK\n")'
```

Expected: `Rd parses OK`

- [ ] **Step 6: Commit**

```bash
git add R/sse-control.R man/runSSEControl.Rd tests/testthat/test-control.R
git commit -m "feat: add covarianceDraw control for joint vs independent_iw OMEGA draws"
```

---

### Task 9b: Implement the weak-identification warning

**Files:**
- Modify: `R/sse-control.R` (validate and store `omegaRseWarn`)
- Modify: `R/sse-omega-draw.R` (`.warnWeakOmega()`)
- Modify: `tests/testthat/test-sse-omega-draw.R`

Task 9 adds the `omegaRseWarn` formal and tests it, but nothing validates the
value, stores it, or ever emits a warning. Without this task the Task 9 test
asserting that a negative value errors cannot pass, and the documented warning
never fires. This task closes that.

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-sse-omega-draw.R`:

```r
test_that("warnWeakOmega fires only above the threshold", {
  omega0 <- diag(c(0.30, 0.30))
  dimnames(omega0) <- list(c("eta.a", "eta.b"), c("eta.a", "eta.b"))

  # RSE 0.20 and 0.20 -- both below 0.5, silent
  expect_silent(
    .warnWeakOmega(omega0, se = c(0.06, 0.06), index = 1:2, threshold = 0.5,
                   mode = "independent_iw")
  )

  # RSE exactly at the threshold -- boundary is NOT "above", so still silent
  expect_silent(
    .warnWeakOmega(omega0, se = c(0.15, 0.06), index = 1:2, threshold = 0.5,
                   mode = "independent_iw")
  )

  # RSE 0.667 on eta.a -- above, warns and NAMES the offending eta
  expect_warning(
    .warnWeakOmega(omega0, se = c(0.20, 0.06), index = 1:2, threshold = 0.5,
                   mode = "independent_iw"),
    "eta\\.a"
  )
})

test_that("warnWeakOmega reports the mode-specific consequence", {
  omega0 <- matrix(0.30, 1L, 1L, dimnames = list("eta.a", "eta.a"))

  # joint, 1x1: exp(SE^2 / (2*omega^2)) is exact here, so claim the inflation
  w <- tryCatch(
    .warnWeakOmega(omega0, se = 0.30, index = 1L, threshold = 0.5,
                   mode = "joint"),
    warning = function(x) conditionMessage(x)
  )
  expect_match(w, "expected mean inflation")

  # independent_iw: report the resulting nu instead
  w2 <- tryCatch(
    .warnWeakOmega(omega0, se = 0.30, index = 1L, threshold = 0.5,
                   mode = "independent_iw"),
    warning = function(x) conditionMessage(x)
  )
  expect_match(w2, "nu")
})

test_that("warnWeakOmega calls the joint figure a proxy for p > 1", {
  # In a correlated block the diagonal is a sum of squared Cholesky elements,
  # so the scalar formula is no longer exact and must not be presented as the
  # actual inflation.
  omega0 <- matrix(c(0.30, 0.05, 0.05, 0.30), 2L, 2L,
                   dimnames = list(c("eta.a", "eta.b"), c("eta.a", "eta.b")))
  w <- tryCatch(
    .warnWeakOmega(omega0, se = c(0.30, 0.06), index = 1:2, threshold = 0.5,
                   mode = "joint"),
    warning = function(x) conditionMessage(x)
  )
  expect_match(w, "proxy")
  expect_no_match(w, "expected mean inflation")
})

test_that("warnWeakOmega ignores unusable standard errors", {
  omega0 <- matrix(0.30, 1L, 1L, dimnames = list("eta.a", "eta.a"))
  expect_silent(
    .warnWeakOmega(omega0, se = NA_real_, index = 1L, threshold = 0.5,
                   mode = "joint")
  )
  expect_silent(
    .warnWeakOmega(omega0, se = 0, index = 1L, threshold = 0.5, mode = "joint")
  )
})
```

- [ ] **Step 2: Run to verify they fail**

Expected: FAIL — `could not find function ".warnWeakOmega"`.

- [ ] **Step 3: Validate and store `omegaRseWarn`**

In `R/sse-control.R`, alongside the other validations, add:

```r
  checkmate::assertNumber(
    omegaRseWarn,
    lower = 0,
    finite = TRUE,
    .var.name = "omegaRseWarn"
  )
```

and add it to the returned `structure()` list, after `covarianceDraw`:

```r
      omegaRseWarn = omegaRseWarn,
```

Gate it to covariance mode alongside `covarianceDraw`:

```r
  if (!omegaRseWarnMissing && parameterSource != "covariance") {
    .abortSSE(
      "{.arg omegaRseWarn} requires {.arg parameterSource = \"covariance\"}."
    )
  }
```

capturing `omegaRseWarnMissing <- missing(omegaRseWarn)` before the value is
touched, exactly as `covarianceDraw` does.

- [ ] **Step 4: Implement the warning**

Append to `R/sse-omega-draw.R`:

```r
#' Warn about weakly identified OMEGA variances in a drawn block
#'
#' `nu = p + 3 + 2*(omega/se)^2` is always greater than `p + 3`, so a threshold
#' on `nu` cannot detect ordinary weak identification. Test the relative
#' standard error directly instead.
#'
#' The consequence differs by mode, so the message says which: for `"joint"`
#' the exact expected mean inflation `exp(se^2 / (2*omega^2))` for a 1x1 block,
#' or a scalar-RSE bias proxy for larger blocks (where a diagonal element is a
#' sum of squared Cholesky elements, so the scalar formula no longer holds);
#' for `"independent_iw"` the
#' resulting degrees of freedom.
#'
#' @param omega0 the block's fitted OMEGA sub-matrix
#' @param se reported standard errors for the block's diagonal
#' @param index the block's eta indices, used only for the message
#' @param threshold `control$omegaRseWarn`
#' @param mode `control$covarianceDraw`
#' @return nothing, called for the warning
#' @noRd
.warnWeakOmega <- function(omega0, se, index, threshold, mode) {
  variances <- diag(omega0)
  etaNames <- rownames(omega0)
  usable <- !is.na(se) & is.finite(se) & se > 0 & variances > 0

  if (!any(usable)) {
    return(invisible(NULL))
  }

  rse <- rep(NA_real_, length(variances))
  rse[usable] <- se[usable] / variances[usable]
  # strictly above: a value exactly at the threshold is not "above" it
  flagged <- which(!is.na(rse) & rse > threshold)

  if (length(flagged) == 0L) {
    return(invisible(NULL))
  }

  detail <- vapply(flagged, function(i) {
    if (identical(mode, "joint")) {
      # exp(RSE^2/2) - 1 is EXACT only for a 1x1 block, where Omega* is
      # lognormal. In a larger block a diagonal element is a sum of squared
      # Cholesky elements, and its induced mean depends on the whole
      # transformed covariance -- so for p > 1 this is a scalar proxy, not the
      # actual inflation. Label it accordingly rather than overstating it.
      inflation <- exp(se[[i]]^2 / (2 * variances[[i]]^2)) - 1
      label <- if (nrow(omega0) > 1L) {
        "scalar-RSE bias proxy"
      } else {
        "expected mean inflation"
      }
      sprintf(
        "%s: relative SE %.0f%%, %s %.0f%%",
        etaNames[[i]], 100 * rse[[i]], label, 100 * inflation
      )
    } else {
      p <- nrow(omega0)
      nu <- p + 3 + 2 * (variances[[i]] / se[[i]])^2
      sprintf(
        "%s: relative SE %.0f%%, implying nu = %.1f",
        etaNames[[i]], 100 * rse[[i]], nu
      )
    }
  }, character(1))

  cli::cli_warn(c(
    "!" = "Weakly identified OMEGA variance{?s} in the block {.val {paste(etaNames, collapse = ', ')}}.",
    "i" = "{detail}",
    "i" = "Parameter uncertainty drawn from these estimates may not be trustworthy."
  ))
  invisible(NULL)
}
```

- [ ] **Step 5: Confirm the call site is present**

`.warnWeakOmega()` must be called once per drawable block, immediately after
`drawn_blocks` is established and before any draw is taken.

**Do not add that loop here.** Task 10 replaces
`.resolveCovarianceParameterSets()` wholesale, so a loop inserted at this point
would be silently deleted when Task 10 lands. The loop is already written into
Task 10's replacement body — verify it is there:

```bash
grep -n "warnWeakOmega" R/sse-helpers.R
```

If Task 10 has already been executed, this returns the call site and there is
nothing to do. If Task 10 has not yet run, this returns nothing — that is
expected, and Task 10 will supply it. Either way, do not hand-insert it.

Warn once per run, not once per replicate — the same block would otherwise
produce `samples` identical warnings.

- [ ] **Step 6: Run to verify the tests pass**

Expected: PASS, including the Task 9 test that a negative `omegaRseWarn`
errors, which only passes once Step 3 lands.

- [ ] **Step 7: Commit**

```bash
git add R/sse-control.R R/sse-omega-draw.R tests/testthat/test-sse-omega-draw.R
git commit -m "feat: warn about weakly identified OMEGA variances"
```

---

### Task 10: Dispatch the per-replicate draw on `covarianceDraw`

**Files:**
- Modify: `R/sse-helpers.R` (`.resolveCovarianceParameterSets()`)
- Modify: `tests/testthat/test-sse-omega-draw.R`

- [ ] **Step 1: Write the failing test**

Append to `tests/testthat/test-sse-omega-draw.R`:

```r
test_that("covariance parameter sets draw both theta and omega", {
  tmp <- tempfile("nlmixr2sse-covdraw-")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE, force = TRUE), add = TRUE)
  nlmixr2utils::withRunSeed(tmp, seed = 11, prefix = "sse")

  nm <- c("tka", "om.eta.ka")
  cov <- diag(c(0.04, 0.001))
  dimnames(cov) <- list(nm, nm)
  omega <- matrix(0.3, 1L, 1L, dimnames = list("eta.ka", "eta.ka"))

  fit <- list(
    theta = c(tka = 0.45),
    omega = omega,
    sigma = matrix(numeric(0), 0L, 0L),
    cov = cov,
    ui = list(iniDf = data.frame(
      name = c("tka", "eta.ka"),
      ntheta = c(1L, NA),
      neta1 = c(NA, 1L),
      neta2 = c(NA, 1L),
      fix = c(FALSE, FALSE),
      stringsAsFactors = FALSE
    ))
  )

  schema <- list(
    thetaCols = "tka",
    omegaCols = "omega(eta.ka,eta.ka)",
    sigmaCols = character(0)
  )

  # NOTE: the loop below runs BOTH modes. Do not collapse it to the default
  # mode only -- a theta-only regression in "joint" is invisible from
  # "independent_iw" and vice versa.

  # both modes must draw theta AND omega
  for (mode in c("joint", "independent_iw")) {
    res <- .resolveCovarianceParameterSets(
      fit = fit,
      samples = 3L,
      outputDir = tmp,
      schema = schema,
      control = runSSEControl(
        parameterSource = "covariance",
        covarianceDraw = mode
      )
    )

    expect_length(res$records, 3L)

    thetas <- vapply(res$records, function(r) unname(r$theta[["tka"]]), numeric(1))
    omegas <- vapply(res$records, function(r) r$omega[1L, 1L], numeric(1))

    # both vary across replicates
    expect_gt(length(unique(thetas)), 1L)
    expect_gt(length(unique(omegas)), 1L)
    # every drawn omega is positive
    expect_true(all(omegas > 0), info = mode)
    # the partition names omega as drawn, not fixed
    expect_true(
      "omega(eta.ka,eta.ka)" %in% res$info$parameterPartition$drawn,
      info = mode
    )
    # the run record states which mode produced the draws
    expect_equal(res$info$covarianceDraw, mode)
  }
})
```

Then append a test that the held-block reason **survives in the returned
object**, not merely on the console:

```r
test_that("held OMEGA blocks and their reasons survive in res$info", {
  # A 2x2 block with a missing off-diagonal: coverage holds the whole block,
  # and the run record must be able to say why after the run has finished.
  tmp <- tempfile("nlmixr2sse-held-")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE, force = TRUE), add = TRUE)
  nlmixr2utils::withRunSeed(tmp, seed = 3, prefix = "sse")

  nm <- c("tka", "om.eta.a", "om.eta.b")   # cov.eta.b.eta.a deliberately absent
  cov <- diag(c(0.04, 0.001, 0.001))
  dimnames(cov) <- list(nm, nm)
  omega <- matrix(c(0.30, 0.05, 0.05, 0.12), 2L, 2L,
                  dimnames = list(c("eta.a", "eta.b"), c("eta.a", "eta.b")))

  fit <- list(
    theta = c(tka = 0.45),
    omega = omega,
    sigma = matrix(numeric(0), 0L, 0L),
    cov = cov,
    ui = list(iniDf = data.frame(
      name = c("tka", "eta.a", "(eta.a,eta.b)", "eta.b"),
      ntheta = c(1L, NA, NA, NA),
      neta1 = c(NA, 1L, 2L, 2L),
      neta2 = c(NA, 1L, 1L, 2L),
      fix = c(FALSE, FALSE, FALSE, FALSE),
      stringsAsFactors = FALSE
    ))
  )

  res <- .resolveCovarianceParameterSets(
    fit = fit,
    samples = 2L,
    outputDir = tmp,
    schema = list(
      thetaCols = "tka",
      omegaCols = c("omega(eta.a,eta.a)", "omega(eta.b,eta.a)",
                    "omega(eta.b,eta.b)"),
      sigmaCols = character(0)
    ),
    control = runSSEControl(parameterSource = "covariance")
  )

  held <- res$info$parameterPartition$heldOmegaBlocks
  expect_length(held, 1L)
  expect_setequal(held[[1L]]$etas, c("eta.a", "eta.b"))
  expect_match(held[[1L]]$reason, "cov\\.eta\\.b\\.eta\\.a")

  # and the block really did stay fitted
  omegas <- vapply(res$records, function(r) r$omega[1L, 1L], numeric(1))
  expect_equal(omegas, rep(0.30, 2L))
})
```

Then append the theta-only regression test — this is the one that catches the
`covFull = FALSE` case, and it must cover **both** modes:

```r
test_that("theta-only fit$cov still draws thetas in BOTH modes", {
  # foceiControl(covFull = FALSE), and the installer's fallback when the full
  # covariance is unavailable, both produce a theta-only fit$cov. Thetas must
  # still be drawn; all OMEGA is held. Guards against gating the whole draw on
  # OMEGA drawability.
  nm <- c("tka", "tcl")
  cov <- diag(c(0.04, 0.09))
  dimnames(cov) <- list(nm, nm)
  omega <- matrix(0.3, 1L, 1L, dimnames = list("eta.ka", "eta.ka"))

  fit <- list(
    theta = c(tka = 0.45, tcl = 1.0),
    omega = omega,
    sigma = matrix(numeric(0), 0L, 0L),
    cov = cov,
    ui = list(iniDf = data.frame(
      name = c("tka", "tcl", "eta.ka"),
      ntheta = c(1L, 2L, NA),
      neta1 = c(NA, NA, 1L),
      neta2 = c(NA, NA, 1L),
      fix = c(FALSE, FALSE, FALSE),
      stringsAsFactors = FALSE
    ))
  )

  schema <- list(
    thetaCols = c("tka", "tcl"),
    omegaCols = "omega(eta.ka,eta.ka)",
    sigmaCols = character(0)
  )

  for (mode in c("independent_iw", "joint")) {
    tmp <- tempfile(paste0("nlmixr2sse-thetaonly-", mode, "-"))
    dir.create(tmp)
    on.exit(unlink(tmp, recursive = TRUE, force = TRUE), add = TRUE)
    nlmixr2utils::withRunSeed(tmp, seed = 5, prefix = "sse")

    res <- .resolveCovarianceParameterSets(
      fit = fit,
      samples = 4L,
      outputDir = tmp,
      schema = schema,
      control = runSSEControl(
        parameterSource = "covariance",
        covarianceDraw = mode
      )
    )

    thetas <- vapply(res$records, function(r) unname(r$theta[["tka"]]), numeric(1))
    omegas <- vapply(res$records, function(r) r$omega[1L, 1L], numeric(1))

    # thetas ARE drawn even with no drawable OMEGA block
    expect_gt(length(unique(signif(thetas, 8))), 1L)
    # OMEGA is held at its fitted value
    expect_equal(omegas, rep(0.3, 4L), info = mode)
    # and reported as fixed, not drawn
    expect_true("tka" %in% res$info$parameterPartition$drawn, info = mode)
    expect_true(
      "omega(eta.ka,eta.ka)" %in% res$info$parameterPartition$fixed,
      info = mode
    )
  }
})
```

- [ ] **Step 2: Run to verify it fails**

Expected: FAIL. The first test fails because every `omega` is currently
identical. The theta-only test fails **for `"joint"` specifically**: with no
drawable OMEGA block the joint spec is skipped, so `thetas` come back constant
and `expect_gt(length(unique(...)), 1L)` fails. It passes for
`"independent_iw"`, which is exactly why the mode loop is required.

- [ ] **Step 3: Implement**

In `R/sse-helpers.R`, replace the whole `.resolveCovarianceParameterSets()` function with:

```r
.resolveCovarianceParameterSets <- function(
  fit,
  samples,
  outputDir,
  schema,
  control
) {
  # No fallback default here: runSSEControl() already resolves covarianceDraw
  # via match.arg(), so a NULL means the control was built by something other
  # than runSSEControl() and we should not silently invent a mode.
  draw_mode <- control$covarianceDraw
  aligned <- .alignedCovariance(fit)
  base_params <- .paramSetFromFit(fit)

  has_theta_draw <- length(aligned$drawNames) > 0L
  chol_cov <- if (has_theta_draw) chol(aligned$cov) else NULL

  omega0 <- base_params$omega
  n_eta <- if (is.matrix(omega0)) nrow(omega0) else 0L
  blocks <- .omegaBlocks(aligned$omegaEntries, n_eta)
  omega_se <- aligned$omegaSe

  # which omega entries actually vary: those in a block with at least one
  # usable standard error
  # Coverage policy: a block is drawn ONLY when every declared entry in it is
  # unfixed and present in fit$cov. Both modes consume this same result -- do
  # NOT re-derive drawability from standard errors here.
  coverage <- .drawableOmegaBlocks(
    blocks,
    aligned$omegaEntries,
    rownames(aligned$fullCov)
  )
  drawn_blocks <- coverage$drawable

  # report exactly the entries of the drawn blocks, so the partition matches
  # what actually varies
  entriesOf <- function(idx) {
    inBlock <- aligned$omegaEntries$row %in% idx &
      aligned$omegaEntries$col %in% idx
    .matrixCoordLabelsForEntries(
      aligned$omegaEntries[inBlock, , drop = FALSE]
    )
  }
  drawn_omega_cols <- unlist(
    lapply(drawn_blocks, entriesOf),
    use.names = FALSE
  ) %||% character(0)

  # surface why any block was held, so the run record explains itself
  for (h in coverage$held) {
    cli::cli_inform(c(
      "i" = "OMEGA block {.val {paste(rownames(omega0)[h$index], collapse = ', ')}} held at its fitted values ({h$reason})."
    ))
  }

  # Warn once per run per drawable block about weakly identified variances --
  # NOT once per replicate, which would emit `samples` identical warnings.
  for (idx in drawn_blocks) {
    .warnWeakOmega(
      omega0[idx, idx, drop = FALSE],
      omega_se[idx],
      idx,
      control$omegaRseWarn,
      draw_mode
    )
  }

  # "joint" needs one spec built once: the transformed covariance over
  # c(theta, stacked omega block vectors), restricted to the drawn entries.
  # Build the joint spec whenever there is ANYTHING to draw -- covered thetas,
  # drawable OMEGA blocks, or both. Gating on OMEGA alone would silently stop
  # drawing thetas for a theta-only fit$cov (e.g. foceiControl(covFull =
  # FALSE)), which the coverage policy explicitly supports.
  # .jointDrawSpec() handles an empty `blocks` list: the Jacobian is NULL, so
  # B collapses to diag(nTheta) and Sigma_T is just the theta covariance.
  joint_spec <- NULL
  if (
    identical(draw_mode, "joint") &&
      (has_theta_draw || length(drawn_blocks) > 0L)
  ) {
    joint_names <- c(
      aligned$drawNames,
      unlist(lapply(drawn_blocks, function(idx) {
        inBlock <- aligned$omegaEntries$row %in% idx &
          aligned$omegaEntries$col %in% idx
        aligned$omegaEntries$covName[inBlock]
      }), use.names = FALSE)
    )
    # Defensive only: .drawableOmegaBlocks() already guarantees every entry of
    # a drawable block is present. Reaching this means the coverage policy and
    # this call disagree, which is an internal inconsistency, not user error.
    missing_joint <- setdiff(joint_names, rownames(aligned$fullCov))
    if (length(missing_joint) > 0L) {
      .abortSSE(
        "Internal error: OMEGA entr{?y/ies} {.val {missing_joint}} passed the coverage check but {?is/are} absent from {.arg fit$cov}."
      )
    }
    joint_spec <- .jointDrawSpec(
      theta = aligned$theta[aligned$drawNames],
      blocks = lapply(drawn_blocks, function(idx) {
        list(omega = omega0[idx, idx, drop = FALSE], index = idx)
      }),
      sigma = aligned$fullCov[joint_names, joint_names, drop = FALSE]
    )
  }

  records <- lapply(seq_len(samples), function(sample_id) {
    drawn <- nlmixr2utils::withRunSeed(
      outputDir,
      key = paste0("covariance-", sample_id),
      prefix = "sse",
      expr = {
        if (identical(draw_mode, "joint")) {
          if (is.null(joint_spec)) {
            # nothing is drawable at all: no covered thetas AND no drawable
            # OMEGA block. Everything keeps its fitted value.
            list(theta = numeric(0), omega = omega0)
          } else {
            # THETA and OMEGA drawn TOGETHER, incorporating their covariance.
            j <- .drawJoint(joint_spec)
            om <- omega0
            for (k in seq_along(drawn_blocks)) {
              idx <- drawn_blocks[[k]]
              om[idx, idx] <- j$omega[[k]]
            }
            list(theta = unname(j$theta), omega = om)
          }
        } else {
          # "independent_iw": THETA and OMEGA drawn INDEPENDENTLY. This shares NWPRI's
          # broad independence factorization, but is NOT NWPRI -- the OMEGA
          # density differs (see the mode comparison in the design doc).
          theta_draw <- if (has_theta_draw) {
            noise <- stats::rnorm(length(aligned$drawNames))
            aligned$theta[aligned$drawNames] + as.numeric(noise %*% chol_cov)
          } else {
            numeric(0)
          }
          omega_draw <- if (n_eta > 0L) {
            # drawn_blocks, NOT blocks: .drawOmega() requires the approved
            # drawable list. Passing the raw block list would redraw blocks the
            # coverage policy held fixed, mutating fixed or uncovered entries.
            .drawOmega(omega0, drawn_blocks, omega_se)
          } else {
            omega0
          }
          list(theta = theta_draw, omega = omega_draw)
        }
      }
    )

    theta <- aligned$theta
    if (has_theta_draw && length(drawn$theta) > 0L) {
      theta[aligned$drawNames] <- drawn$theta
    }

    list(
      sample = as.integer(sample_id),
      sourceSample = NA_integer_,
      source = "covariance",
      hypothesis = "simulation",
      modelLabel = "simulation",
      role = "simulation",
      theta = theta,
      omega = drawn$omega,
      sigma = base_params$sigma
    )
  })

  list(
    records = records,
    info = list(
      mode = "covariance",
      covarianceDraw = draw_mode,
      parameterPartition = list(
        drawn = unname(c(aligned$drawNames, drawn_omega_cols)),
        fixed = unname(c(
          setdiff(schema$thetaCols, aligned$drawNames),
          setdiff(schema$omegaCols, drawn_omega_cols),
          schema$sigmaCols
        )),
        # Persist WHY each block was held. cli_inform() above only reaches the
        # console; without this the reason is lost the moment the run ends, and
        # a saved run cannot explain why its OMEGA did not vary.
        heldOmegaBlocks = lapply(coverage$held, function(h) {
          list(
            etas = rownames(omega0)[h$index],
            reason = h$reason
          )
        })
      ),
      estimationInitialValues = "reference_fit"
    )
  )
}
```

This references a small helper that converts entry-table rows into the canonical `omega(<eta>,<eta>)` raw-results column labels. Add it to `R/sse-omega-draw.R`:

```r
#' Canonical raw-results column labels for OMEGA entries
#'
#' The raw-results schema labels OMEGA entries `omega(<rowEta>,<colEta>)`,
#' which is a different convention from `fit$cov`'s `om.`/`cov.` naming.
#' @noRd
.matrixCoordLabelsForEntries <- function(entries) {
  if (nrow(entries) == 0L) {
    return(character(0))
  }
  paste0("omega(", entries$rowName, ",", entries$colName, ")")
}
```

- [ ] **Step 4: Thread `control` through the caller**

`.resolveCovarianceParameterSets()` now takes `control`, so its caller must
pass it. In `R/sse-helpers.R`, find the call inside
`.resolveSimulationParameters()`:

```r
  .resolveCovarianceParameterSets(
    fit = fit,
    samples = samples,
    outputDir = outputDir,
    schema = schema
  )
```

Replace with:

```r
  .resolveCovarianceParameterSets(
    fit = fit,
    samples = samples,
    outputDir = outputDir,
    schema = schema,
    control = control
  )
```

`.resolveSimulationParameters()` already receives `control`, so nothing further
up the chain changes.

- [ ] **Step 5: Run to verify it passes**

Expected: PASS.

- [ ] **Step 6: Verify end to end against a real fit**

This is the actual user-facing defect. Run:

```bash
NOT_CRAN=true Rscript -e '
devtools::load_all(quiet = TRUE)
suppressMessages(library(nlmixr2))
ref <- function() {
  ini({
    tka <- log(1.57); tcl <- log(2.72); tv <- log(31.5)
    eta.ka ~ 0.6
    add.sd <- 0.7
  })
  model({
    ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv)
    cp <- linCmt(); cp ~ add(add.sd)
  })
}
fit <- suppressMessages(suppressWarnings(nlmixr2(ref, nlmixr2data::theo_sd,
  est = "focei", control = list(print = 0L, covMethod = "r"))))
tmp <- tempfile("sse-cov-e2e-")
res <- suppressMessages(runSSE(fit, samples = 3L, seed = 42,
  control = runSSEControl(workers = 1L, parameterSource = "covariance"),
  outputDir = tmp, restart = TRUE))
cat("class:", paste(class(res), collapse = ","), "\n")
cat("drawn:", paste(res$runInfo$parameterSourceInfo$parameterPartition$drawn,
                    collapse = ", "), "\n")
'
```

Expected: completes without the `contains theta name "om.eta.ka"` error, and the drawn list includes an `omega(...)` entry alongside the thetas.

Then repeat the same call with `covarianceDraw = "independent_iw"` added to the
`runSSEControl(...)` and confirm it also completes.

- [ ] **Step 7: Run the full suite**

Expected: FAIL still 5, PASS higher.

- [ ] **Step 8: Commit**

```bash
git add R/sse-helpers.R R/sse-omega-draw.R tests/testthat/test-sse-omega-draw.R
git commit -m "feat: draw OMEGA with uncertainty in covariance parameter mode"
```

---

### Task 11: Document

**Files:**
- Modify: `NEWS.md`
- Modify: `README.md`

- [ ] **Step 1: Add the NEWS entry**

At the top of `NEWS.md`, immediately under the `# nlmixr2sse 0.1` heading, add:

```markdown
* `parameterSource = "covariance"` now draws OMEGA as well as THETA. Recent
  `nlmixr2est` reports a joint fixed/random-effect covariance in `fit$cov`
  (OMEGA entries named `om.<eta>` / `cov.<eta>.<eta>`), which previously made
  this mode abort with "contains theta name ... not present in `fit$theta`".
* New `runSSEControl(covarianceDraw = )` selects how the draw is taken:
  * `"independent_iw"` (default) draws THETA multivariate-Normal and OMEGA from
    a mean-centred inverse-Wishart per OMEGA block, independently of each
    other, with degrees of freedom moment-matched to the reported OMEGA
    standard errors. The OMEGA draw is delegated to `rxode2::cvPost()`.
  * `"joint"` instead draws THETA and OMEGA **together** from `fit$cov` on a
    log-Cholesky-transformed scale, incorporating the estimated THETA/OMEGA
    covariance through a first-order delta approximation.

  Both give positive-definite OMEGA draws. **Neither implements NONMEM's
  `$PRIOR NWPRI`, and neither claims PsN parity.**
* `"joint"` uses covariance information `"independent_iw"` discards, but its
  non-linear back-transform inflates OMEGA means by roughly
  `exp(SE^2 / (2 * Omega^2))` — about 2% at 20% relative standard error, but
  65% at 100% — so it is opt-in rather than the default. New
  `runSSEControl(omegaRseWarn = )` (default `0.5`) warns when a drawn OMEGA
  variance's relative standard error exceeds the threshold.
* OMEGA blocks are drawn only when `fit$cov` covers every declared element of
  the block. A block containing a fixed or uncovered element is held entirely
  at its fitted values and reported in the run's parameter partition. A
  theta-only `fit$cov` (for example from `foceiControl(covFull = FALSE)`) is
  supported: thetas are drawn and all OMEGA is held fixed.
* **Reproducibility note:** because OMEGA now varies between replicates,
  `parameterSource = "covariance"` runs will not reproduce results from earlier
  versions at the same seed. The two `covarianceDraw` modes also differ from
  each other at the same seed, so the mode is recorded in `run_info` and
  checked when a run is resumed.
* Residual-error uncertainty was already included and is unchanged: in nlmixr2
  residual parameters are thetas, so they are covered by the THETA draw.
  THETA draws remain unconstrained, so a bounded theta or a residual-error
  standard deviation can still be drawn outside its valid range.
```

- [ ] **Step 2: Update the README's parameter-source table**

In `README.md`, find the `covariance` row of the "Three parameter sources" table. It currently reads:

```
| `covariance` | `runSSE(fit, samples = 50, seed = 42, control = runSSEControl(parameterSource = "covariance"))` | Draws theta values from `fit$cov`; OMEGA and SIGMA stay fixed. |
```

Replace the description cell so the row reads:

```
| `covariance` | `runSSE(fit, samples = 50, seed = 42, control = runSSEControl(parameterSource = "covariance"))` | Draws THETA from `fit$cov` and OMEGA from a mean-centred inverse-Wishart per block. Residual error is drawn with the thetas. Set `covarianceDraw = "joint"` to draw THETA and OMEGA together instead. |
```

Then add a short subsection immediately after that table:

```markdown
### Covariance draw modes

`parameterSource = "covariance"` supports two draws, selected with
`runSSEControl(covarianceDraw = )`:

- **`"independent_iw"` (default)** — THETA is drawn multivariate-Normal from
  its own sub-block of `fit$cov`; OMEGA is drawn from a mean-centred
  inverse-Wishart per block, independently. Simple, and free of transformation
  bias.
- **`"joint"`** — THETA and OMEGA are drawn together from `fit$cov`, so the
  covariance between them is incorporated. OMEGA is transformed to a
  log-Cholesky scale for the draw, which makes every sampled matrix
  positive-definite.

Neither mode implements NONMEM's `$PRIOR NWPRI`.

`"joint"` uses information the default discards: the THETA/OMEGA covariance is
small in rich designs but substantial in the sparse ones SSE is usually used to
plan. It is nevertheless opt-in, because its non-linear back-transform inflates
OMEGA means by roughly `exp(SE^2 / (2 * Omega^2))` — about 2% at 20% relative
standard error, but 65% at 100%. Poorly identified variance components sit at
the high end of that range, so `"joint"` warns when a drawn OMEGA variance
exceeds `omegaRseWarn` (default 50% relative standard error).

An OMEGA block is drawn only when `fit$cov` covers all of its declared
elements; otherwise the whole block stays at its fitted values. With a
theta-only `fit$cov`, thetas are still drawn and all OMEGA is held fixed.

THETA draws are unconstrained in both modes, so a bounded theta or a
residual-error standard deviation can be drawn outside its valid range.
```

- [ ] **Step 3: Verify the README claim is accurate**

Confirm no other line in `README.md` still claims OMEGA stays fixed:

```bash
grep -n -i "omega" README.md
```

Expected: no remaining text asserting OMEGA is held fixed in covariance mode. If any is found, correct it.

- [ ] **Step 4: Run the full suite**

Expected: FAIL still 5, PASS unchanged from Task 6.

- [ ] **Step 5: Commit**

```bash
git add NEWS.md README.md
git commit -m "docs: document OMEGA uncertainty in covariance parameter mode"
```

---

### Task 11b: The required edge-case tests

**Files:**
- Modify: `tests/testthat/test-sse-omega-draw.R`
- Modify: `tests/testthat/test-sse-omega-joint.R`

These are not optional robustness extras. Each one guards a specific failure
mode identified in review, and several of them describe bugs that were present
in an earlier draft of this plan. Do not skip any as "unlikely".

Tasks 2 and 10 already cover: a correlated block with a fixed element, a
missing off-diagonal, exclusion-reason reporting, theta-only in both modes, and
the OMEGA-only zero-theta case. This task adds the rest.

- [ ] **Step 1: Write the scale-range and conditioning tests**

Append to `tests/testthat/test-sse-omega-joint.R`:

```r
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
```

- [ ] **Step 2: Write the multi-block ordering test**

Append to `tests/testthat/test-sse-omega-joint.R`:

```r
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
```

- [ ] **Step 3: Run to verify they fail or pass honestly**

Run both files:

```bash
NOT_CRAN=true Rscript -e "devtools::load_all(quiet=TRUE); testthat::test_file('tests/testthat/test-sse-omega-joint.R', reporter='summary')"
NOT_CRAN=true Rscript -e "devtools::load_all(quiet=TRUE); testthat::test_file('tests/testthat/test-sse-omega-draw.R', reporter='summary')"
```

Expected: PASS if Tasks 6-8 were implemented correctly. If the scale-range or
ill-conditioned tests fail, the Jacobian step is not genuinely relative — fix
`.numericJacobian()`, do not relax the test.

- [ ] **Step 4: Commit**

```bash
git add tests/testthat/test-sse-omega-draw.R tests/testthat/test-sse-omega-joint.R
git commit -m "test: cover OMEGA scale range, conditioning, and multi-block ordering"
```

---

### Task 12: Correct the statements this work makes false

**Files:**
- Modify: `R/run-sse.R` (the covariance-mode cli message)
- Modify: `README.md`
- Modify: `vignettes/runSSE.Rmd`
- Modify: `tests/testthat/test-run-sse.R`

Task 11 documents the new behaviour. This task removes the old claims that
contradict it. They are separate tasks because these edits touch a vignette and
a masked test, and reviewing them together with the NEWS entry obscures both.

- [ ] **Step 1: Fix the run banner**

In `R/run-sse.R`, find:

```r
      "i" = "Thetas are drawn from {.field fit$cov}; OMEGA and SIGMA stay at the fitted point estimates."
```

Replace with a message that reports the actual mode and coverage:

```r
      "i" = "Parameters are drawn from {.field fit$cov} ({.field {control$covarianceDraw}} draw).",
      "i" = "OMEGA blocks without full covariance coverage stay at their fitted values."
```

- [ ] **Step 2: Fix the README parameter-source table**

Task 11 already rewrote this row. Verify no other line still claims OMEGA is
held fixed:

```bash
grep -n -i "omega" README.md
```

Expected: no remaining text asserting OMEGA stays fixed in covariance mode.

- [ ] **Step 3: Fix the vignette**

In `vignettes/runSSE.Rmd`, find the covariance-mode passage near "theta-only
uncertainty: OMEGA and SIGMA" and the follow-on "For full uncertainty in OMEGA
and SIGMA as well, use". Rewrite them to state that covariance mode now draws
OMEGA too, that `covarianceDraw` selects how, and that residual error is drawn
with the thetas because nlmixr2 parameterizes it as thetas.

Verify the vignette still renders:

```bash
Rscript -e "rmarkdown::render('vignettes/runSSE.Rmd', quiet = TRUE, output_dir = tempdir())"
```

Expected: exit status 0.

- [ ] **Step 4: Correct the masked assertion**

In `tests/testthat/test-run-sse.R`, find:

```r
  expect_equal(length(unique(signif(cov_omega, 8))), 1L)
```

Replace with:

```r
  expect_gte(length(unique(signif(cov_omega, 8))), 2L)
```

This test currently `skip()`s, so it will not run — correct it anyway so it is
right when the sibling-package defect is fixed.

- [ ] **Step 5: Confirm the theta-only test still passes untouched**

`test-run-sse.R:168` ("covariance parameter resolution draws theta and keeps
omega fixed") must still pass **without modification** — its fixture has a
theta-only `cov`. If it fails, the coverage policy is wrong; fix the policy,
not the test.

Run: `NOT_CRAN=true Rscript -e "devtools::load_all(quiet=TRUE); testthat::test_file('tests/testthat/test-run-sse.R', reporter='summary')"`

- [ ] **Step 6: Commit**

```bash
git add R/run-sse.R README.md vignettes/runSSE.Rmd tests/testthat/test-run-sse.R
git commit -m "docs: correct statements that OMEGA stays fixed in covariance mode"
```

---

### Task 13: Put `covarianceDraw` in the reproducibility contract

**Files:**
- Modify: `R/run-sse.R`
- Modify: `R/recompute-sse.R`
- Modify: `tests/testthat/test-run-sse.R`

Switching modes changes every simulated parameter set at the same seed, so the
mode must be recorded and checked — exactly as `rxThreads` already is.

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-run-sse.R`:

```r
test_that("validateResumeRequest aborts on a covarianceDraw mismatch", {
  run_info <- list(
    fitName = "fake_sse_fit",
    samples = 2L,
    parameterSource = "covariance",
    estimateSimulation = TRUE,
    covarianceDraw = "joint"
  )

  err <- capture_sse_error(
    .validateResumeRequest(
      existingRunInfo = run_info,
      fitName = "fake_sse_fit",
      samples = 2L,
      control = runSSEControl(
        parameterSource = "covariance",
        covarianceDraw = "independent_iw"
      ),
      requestedLabels = character(0)
    )
  )

  expect_s3_class(err, "error")
  expect_match(conditionMessage(err), "covarianceDraw")
  expect_match(conditionMessage(err), "joint")
})

test_that("validateResumeRequest aborts on a legacy covariance run", {
  # No recorded covarianceDraw means the run predates OMEGA drawing, so its
  # replicates hold OMEGA fixed. Resuming would mix two simulation
  # distributions -- this must abort, not warn.
  run_info <- list(
    fitName = "fake_sse_fit",
    samples = 2L,
    parameterSource = "covariance",
    estimateSimulation = TRUE
  )

  err <- capture_sse_error(
    .validateResumeRequest(
      existingRunInfo = run_info,
      fitName = "fake_sse_fit",
      samples = 2L,
      control = runSSEControl(parameterSource = "covariance"),
      requestedLabels = character(0)
    )
  )

  expect_s3_class(err, "error")
  expect_match(conditionMessage(err), "restart")
})

test_that("a recorded NA covarianceDraw is treated as legacy, not a mismatch", {
  # An addModels run against a legacy directory records NA by design. Resuming
  # that directory must give the legacy abort, not a confusing message saying
  # the original mode was "NA".
  run_info <- list(
    fitName = "fake_sse_fit",
    samples = 2L,
    parameterSource = "covariance",
    estimateSimulation = TRUE,
    covarianceDraw = NA_character_
  )

  err <- capture_sse_error(
    .validateResumeRequest(
      existingRunInfo = run_info,
      fitName = "fake_sse_fit",
      samples = 2L,
      control = runSSEControl(parameterSource = "covariance"),
      requestedLabels = character(0)
    )
  )

  expect_s3_class(err, "error")
  expect_match(conditionMessage(err), "restart")
  expect_no_match(conditionMessage(err), "NA")
})

test_that("a legacy run in a non-covariance mode resumes unaffected", {
  # The legacy check is gated on parameterSource == "covariance"; a fixed or
  # rawres run has no OMEGA-draw provenance to be incompatible about.
  run_info <- list(
    fitName = "fake_sse_fit",
    samples = 2L,
    parameterSource = "fixed",
    estimateSimulation = TRUE
  )

  expect_silent(
    .validateResumeRequest(
      existingRunInfo = run_info,
      fitName = "fake_sse_fit",
      samples = 2L,
      control = runSSEControl(parameterSource = "fixed"),
      requestedLabels = character(0)
    )
  )
})

test_that("recordedCovarianceDraw keeps a legacy addModels run as NA", {
  # Exercises the real helper, not a copy of its expression -- a test that
  # re-evaluates `%||%` inline would pass before the implementation existed and
  # could never catch a provenance regression.
  legacy <- list(covarianceDraw = NULL)
  known <- list(covarianceDraw = "joint")
  addCtl <- runSSEControl(
    parameterSource = "covariance",
    addModels = TRUE
  )

  expect_true(is.na(.recordedCovarianceDraw(legacy, addCtl)))
  expect_equal(.recordedCovarianceDraw(known, addCtl), "joint")
})

test_that("recordedCovarianceDraw records NA outside covariance mode", {
  # runSSEControl() always resolves a default, so without this a fixed or
  # rawres run would be stamped "independent_iw" despite never drawing one.
  expect_true(is.na(
    .recordedCovarianceDraw(list(), runSSEControl(parameterSource = "fixed"))
  ))
  expect_equal(
    .recordedCovarianceDraw(
      list(),
      runSSEControl(parameterSource = "covariance")
    ),
    "independent_iw"
  )
})
```

- [ ] **Step 2: Run to verify they fail**

Expected: all of the new tests fail, but for different reasons -- the mismatch and legacy-abort tests because no error is raised yet, and the `.recordedCovarianceDraw()` tests because that helper does not exist. Check each failure message matches its cause; a test failing for an unexpected reason means something else is wrong.

- [ ] **Step 3: Record the mode**

In `R/run-sse.R`, alongside the other `run_info` provenance fields, add:

```r
  run_info$covarianceDraw <- .recordedCovarianceDraw(
    existing_run_info,
    control
  )
```

with the rule extracted into a testable helper in `R/sse-omega-draw.R`:

```r
#' Decide what `covarianceDraw` a run should record
#'
#' Two rules, both about not writing something untrue:
#'
#' 1. An add-models run must not overwrite the value describing how the SAVED
#'    datasets were produced. If the original recorded nothing, the honest
#'    record is `NA` (legacy) -- writing the current mode would claim those
#'    datasets were generated a way they were not.
#' 2. A run that is not in covariance mode never took a covariance draw, so
#'    recording a mode would be meaningless. `runSSEControl()` always resolves
#'    a default, so without this a `fixed` or `rawres` run would be stamped
#'    `"independent_iw"`.
#' @noRd
.recordedCovarianceDraw <- function(existingRunInfo, control) {
  if (isTRUE(control$addModels)) {
    return(existingRunInfo$covarianceDraw %||% NA_character_)
  }
  if (!identical(control$parameterSource, "covariance")) {
    return(NA_character_)
  }
  control$covarianceDraw
}
```

This mirrors the `addModels` preservation pattern used by `parameterSource`,
`estimateSimulation`, and `rxThreads`: an add-models run must not overwrite the
value that describes how the saved datasets were produced. The `NA_character_`
fallback is the difference from `rxThreads` — a missing `covarianceDraw` is not
merely unknown, it identifies a run whose OMEGA was held fixed, and recording
the current mode would misdescribe it.

- [ ] **Step 4: Validate on resume**

In `R/recompute-sse.R`, inside `.validateResumeRequest()`, add alongside the
existing `parameterSource` check:

```r
  if (identical(control$parameterSource, "covariance")) {
    recordedDraw <- existingRunInfo$covarianceDraw
    # NULL and NA both mean "legacy, mode unknown". NA specifically arises
    # after an addModels run against a legacy directory, which records NA by
    # design -- without treating it as legacy here, that run would fall into
    # the generic mismatch branch below and report the original mode as "NA".
    if (is.null(recordedDraw) || is.na(recordedDraw)) {
      # A covariance run with no recorded mode necessarily predates OMEGA
      # drawing, so its existing replicates hold OMEGA at the fitted values.
      # Appending OMEGA-varying replicates to those would put two different
      # simulation distributions in one study. This is NOT the rxThreads case,
      # where the thread count is merely unknown but the draw is the same kind
      # -- here the old behaviour is known and known to be incompatible. Abort.
      .abortSSE(
        paste0(
          "This run directory predates {.arg covarianceDraw} tracking, so its ",
          "existing replicates were simulated with OMEGA held at the fitted ",
          "values. Resuming would append replicates that vary OMEGA, mixing two ",
          "simulation distributions in one study. Use {.code restart = TRUE} to ",
          "start a fresh run."
        )
      )
    } else if (!identical(recordedDraw, control$covarianceDraw)) {
      .abortSSE(
        paste0(
          "Existing run directory used {.arg covarianceDraw = {recordedDraw}}, ",
          "but this run requests {.arg covarianceDraw = {control$covarianceDraw}}. ",
          "The draw mode changes every simulated parameter set, so the replicates would not be comparable. ",
          "Use {.code runSSEControl(covarianceDraw = \"{recordedDraw}\")} to resume, or {.code restart = TRUE} for a fresh run."
        )
      )
    }
  }
```

Gate on `parameterSource == "covariance"` so runs in other modes are
unaffected.

- [ ] **Step 5: Also check it on the `addModels` path**

The `addModels` branch of `runSSE()` does **not** call
`.validateResumeRequest()`. Add a check there alongside the existing
`updateFix` validation, warning (not aborting) on a mismatch, since add-models
refits against saved datasets:

```r
    if (identical(control$parameterSource, "covariance")) {
      recordedDraw <- existing_run_info$covarianceDraw
      if (is.null(recordedDraw) || is.na(recordedDraw)) {
        # Legacy directory: the mode is genuinely unknown, so do not phrase
        # this as a mismatch against "NA". Same NULL-or-NA legacy test as the
        # resume path, but a warning rather than an abort, because add-models
        # simulates nothing new.
        cli::cli_warn(c(
          "!" = "The original run predates {.arg covarianceDraw} tracking, so the draw mode behind its saved datasets is unknown.",
          "i" = "Those datasets are reused unchanged; the recorded {.field covarianceDraw} stays unset."
        ))
      } else if (!identical(recordedDraw, control$covarianceDraw)) {
        cli::cli_warn(c(
          "!" = "The original run used {.arg covarianceDraw = {recordedDraw}}, but this add-models run requests {.arg covarianceDraw = {control$covarianceDraw}}.",
          "i" = "Saved simulated datasets are reused unchanged; the recorded {.field covarianceDraw} keeps the original value."
        ))
      }
    }
```

- [ ] **Step 6: Run to verify the tests pass**

Expected: PASS.

- [ ] **Step 7: Run the full suite**

Expected: no new failures beyond the documented list.

- [ ] **Step 8: Commit**

```bash
git add R/run-sse.R R/recompute-sse.R tests/testthat/test-run-sse.R
git commit -m "feat: record and validate covarianceDraw across resume and addModels"
```

---

## Self-Review

**Task map:** 1-2 shared plumbing · 3-4 `"independent_iw"` · 5 partition (shared, fixes
the defect) · 6-8 `"joint"` · 9 control argument · 10 dispatch · 11 docs.

**Spec coverage:**

- `om.`/`cov.` name construction generated locally and intersected, not parsed — Task 1.
- Per-block structure via connected components — Task 2 (Steps 1-5).
- **Coverage policy** — `.drawableOmegaBlocks()` in Task 2 (Steps 6-10). A block is drawn only when every declared entry is unfixed and present in `fit$cov`; otherwise the whole block is held and the reason recorded. Both modes consume this single result (Task 10), neither re-derives drawability from standard errors.
- **Theta-only `fit$cov` draws thetas in BOTH modes** — Task 10. The joint spec is built whenever there are covered thetas *or* drawable blocks, never gated on OMEGA alone, and the regression test loops over both modes because a `"joint"`-only failure is invisible from the default mode.
- `"independent_iw"`: OMEGA inverse-Wishart centred on `fit$omega`; `nu` moment-matched per block via `min(nu_i)` — Tasks 3, 4.
- Defect fixed (no abort on OMEGA names) — Task 5.
- `"joint"`: log-Cholesky transform, PD by construction — Task 6.
- `"joint"`: delta-method Jacobian carrying `fit$cov` onto the transformed scale — Task 7.
- `"joint"`: `Σ_T = B Σ Bᵀ`, single joint draw, cross-terms incorporated to first order — Task 8.
- `covarianceDraw` control, defaulting to `"independent_iw"`, gated to `parameterSource = "covariance"` — Task 9.
- THETA Normal draw from the theta sub-block (`"independent_iw"` path) — Task 5 + Task 10.
- Fixed/uncovered parameters held at fitted values — Task 4 (`NULL` spec skipped), Task 5 (fixed diagonals skipped), Task 10 (`drawn_blocks` filter).
- `parameterPartition` extended to OMEGA; `covarianceDraw` recorded in run info — Task 10.
- No SIGMA draw; residual error rides the THETA draw — stated in Background, asserted in Task 11 docs.
- Error handling: unusable `nu`, degenerate variance, singular `Psi` (Tasks 3, 4); unknown `fit$cov` name (Task 5); non-PD `Σ_T` and missing joint entries (Tasks 8, 10); `covarianceDraw` misuse (Task 9).
- Statistical assertions: `"independent_iw"` mean/SD recovery (Task 4); `"joint"` correlation recovery and universal PD (Task 8).
- Reproducibility break and the joint-vs-independent_iw trade-off documented — Task 11.

**Placeholder scan:** No TBD/TODO. Every code step shows complete code; every run step shows the command and expected output.

**Type consistency:** `.omegaEntryTable()` returns columns `row/col/rowName/colName/fix/covName/diagonal`, used consistently by `.omegaBlocks()` (`row`, `col`, `diagonal`), Task 5 (`covName`, `diagonal`, `rowName`, `fix`), and Task 6 (`row`, `col`). `.omegaWishartSpec()` returns `list(omega0, nu, p)` -- the fitted block, NOT a precomputed Psi, because `cvPost()` pre-scaling happens at draw time -- consumed by `.drawOmegaBlock()` (`omega0`, `nu`, `p`) and Task 3's test (`omega0`, `nu`, `p`). `.alignedCovariance()` returns `theta/cov/fullCov/drawNames/omegaSe/omegaEntries`; `cov` is the theta-only sub-block while `fullCov` is the whole matrix that `"joint"` needs for its cross-terms, and `drawNames` keeps its original meaning (theta names only).

**Known risk carried from the spec:** the `om.`/`cov.` naming convention comes from unexported `nlmixr2est` internals. Task 1 Step 5 verifies it against a live fit and instructs the implementer to STOP if it no longer holds, rather than silently misparsing.
