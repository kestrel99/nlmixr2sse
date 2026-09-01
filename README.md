# nlmixr2sse

`nlmixr2sse` implements stochastic simulation and estimation (SSE) using
`nlmixr2` for parameter-recovery exercises, model-comparison sample-size work,
and Type I / power evaluation.

## What it does

- Generate `N` simulated datasets from a fitted `nlmixr2` reference model.
- Refit the simulation model and any named alternatives to each replicate.
- Summarize parameter bias, RMSE, and OFV-based power signals from canonical
  raw-results output.

## Relationship to sibling packages

`nlmixr2sse` reuses the canonical raw-results schema and shared run-cache
helpers from `nlmixr2utils`, so `raw_results.*` files move directly between
`nlmixr2sse`, `nlmixr2boot`, and `nlmixr2sir` (and any other packages that
rely on `nlmixr2utils`) without format conversion.

## Install

```r
remotes::install_github("kestrel99/nlmixr2sse")
```

## Vignettes

For detailed definitions, assumptions, statistical targets, and implementation
notes, see the [technical reference](docs/sse-technical-reference.md).

Running SSE end to end: restart and resume, plotting, and
uncertainty-driven parameter sources.

```r
vignette("runSSE", package = "nlmixr2sse")
```

A worked sample-size evaluation: how many subjects are needed to detect
a covariate effect?

```r
vignette("sse-power", package = "nlmixr2sse")
```

## End-to-end example

```r
library(nlmixr2)
library(nlmixr2data)
library(nlmixr2sse)

ref_model <- function() {
  ini({
    tka <- log(1.57)
    tcl <- log(2.72)
    tv <- log(31.5)
    eta.ka ~ 0.6
    add.sd <- 0.7
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl)
    v <- exp(tv)
    cp <- linCmt()
    cp ~ add(add.sd)
  })
}

alt_model <- function() {
  ini({
    tka <- log(1.57)
    tcl <- log(2.72)
    tv <- log(31.5)
    add.sd <- 0.7
  })
  model({
    ka <- exp(tka)
    cl <- exp(tcl)
    v <- exp(tv)
    cp <- linCmt()
    cp ~ add(add.sd)
  })
}

fit <- nlmixr2(
  ref_model,
  nlmixr2data::theo_sd,
  est = "focei",
  control = list(print = 0L, covMethod = "r")
)
alt_fit <- nlmixr2(
  alt_model,
  nlmixr2data::theo_sd,
  est = "focei",
  control = list(print = 0L, covMethod = "r")
)

comparison <- sseComparison(
  full = "simulation",
  reduced = "no_eta",
  criticalValue = stats::qchisq(0.90, df = 1),
  label = "eta on ka"
)

sse <- runSSE(
  fit,
  alternativeModels = sseModel(alt_fit, label = "no_eta"),
  comparisons = comparison,
  samples = 50,
  seed = 42
)

summary(sse)
comparisonSummary(sse)
```

This example tests a single variance component at the boundary of its parameter
space. Its illustrative cutoff is the 5% cutoff from the commonly used 50:50
mixture of point mass at zero and chi-square with one degree of freedom. Confirm
that a boundary reference is appropriate for the model and estimation method
before using it in an analysis. Because this comparison uses a custom critical
value rather than an ordinary chi-square `df`, it is not eligible for
distribution-based PPE.

Use `samples = 50` only for a quick example. Choose the number of replicates
from the Monte Carlo precision the analysis needs, rather than from a fixed
rule of thumb. For an estimated probability `p`, the approximate Monte Carlo
standard error is `sqrt(p * (1 - p) / n)`. At `p = 0.05`, for example, it is
about 2.2 percentage points with 100 paired-evaluable replicates and 0.7
percentage points with 1000. `comparisonSummary()` reports both the effective
paired denominator and an exact binomial confidence interval; inspect that
interval especially when the observed rate is zero or one.

## Explicit model comparisons

Define important hypotheses explicitly with `sseComparison()`. The test
statistic is always

```text
T = OFV(reduced) - OFV(full)
```

and the null is rejected when `T` exceeds the comparison's critical value.
Using `"simulation"` names the data-generating model regardless of the label
derived from `fit`. Supplying `df` lets `sseComparison()` calculate the
chi-square cutoff from `alpha`; supplying `criticalValue` instead supports a
custom reference but disables distribution-based PPE.

If `comparisons` is omitted, the package can construct legacy
simulation-versus-alternative comparisons, but it must infer degrees of freedom
from raw-results parameter columns. That count can be wrong for fixed,
constrained, or non-nested hypotheses, so explicit comparisons are recommended.

`comparisonSummary()` calculates the rejection rate only among replicates in
which both models have finite, accepted OFVs. It reports attempted, evaluable,
and excluded counts, the Monte Carlo standard error, and an exact binomial
interval. If fit failure or filtering is related to the simulated data, this
complete-case rate is conditional and can differ from the unconditional design
operating characteristic.

## Three parameter sources

SSE parameters can be obtained in three ways.

| Mode | One-line example | Notes |
| --- | --- | --- |
| `fixed` | `runSSE(fit, samples = 50, seed = 42)` | Reuses the fitted point estimates for every replicate. |
| `covariance` | `runSSE(fit, samples = 50, seed = 42, control = runSSEControl(parameterSource = "covariance"))` | Draws THETA from `fit$cov` and OMEGA from a mean-centred inverse-Wishart per block. Residual-error parameters represented as THETAs are drawn with them. Set `covarianceDraw = "joint"` to draw THETA and OMEGA together instead. |
| `rawres` | `runSSE(fit, samples = 50, seed = 42, control = runSSEControl(parameterSource = "rawres", rawresInput = file.path(boot_dir, "raw_results.csv")))` | Uses one canonical raw-results row per replicate, typically from `nlmixr2boot::runBootstrap()`. |

### Covariance draw modes

`parameterSource = "covariance"` supports two draws, selected with
`runSSEControl(covarianceDraw = )`:

- **`"independent_iw"` (default)** — THETA is drawn multivariate-Normal from
  its own sub-block of `fit$cov`; OMEGA is drawn from a mean-centred
  inverse-Wishart per block, independently. Simple, and free of transformation
  bias. A single inverse-Wishart degrees of freedom is chosen per block, so
  only its binding diagonal element matches its reported standard error
  exactly; other elements can be over-dispersed.
- **`"joint"`** — THETA and OMEGA are drawn together from `fit$cov`, so the
  covariance between them is incorporated through a first-order delta
  approximation. OMEGA is transformed to a log-Cholesky scale for the draw,
  which makes every sampled matrix positive-definite.

For those familiar with PsN's approach to SSE, neither of these approaches
implements anything similar to NONMEM's `$PRIOR NWPRI`.

`"joint"` uses information that the default would normally discard, the full
THETA/OMEGA covariance. This is usually small in rich designs but can be
substantial in the sparse ones SSE is often used to plan. It is nevertheless
not the default, because its non-linear back-transform inflates
OMEGA means by roughly `exp(SE^2 / (2 * Omega^2))` — about 2% at 20% relative
standard error, but 65% at 100%. Poorly identified variance components
typically occupy the high end of that range. A drawable OMEGA block triggers
one weak-identification warning when a fitted diagonal element has reported
`SE / OMEGA` greater than `omegaRseWarn` (default 0.5). The check uses the
fitted value and its reported standard error; it is not a limit applied to
individual draws.

An OMEGA block is drawn only when `fit$cov` covers all of its declared
elements; otherwise the whole block stays at its fitted values. With a
THETA-only `fit$cov`, THETAs are still drawn but OMEGA is fixed.

THETA draws are unconstrained in both modes, so a bounded THETA or a
residual-error standard deviation can be drawn outside its valid range.

For a completed covariance-mode run, inspect the values actually drawn rather
than assuming the approximation reproduced every target moment:

```r
parameterDrawSummary(covariance_sse)
plotSSEParameterDraws(covariance_sse)
```

The summary reports realized versus target moments, out-of-domain THETA draws,
OMEGA positive-definiteness checks, joint-mode raw-scale mean drift, and
independent-inverse-Wishart dispersion ratios.

## Resume and recompute

The core output from `runSSE()` includes `run_info.rds`, `sse_state.rds`,
`initial_estimates.csv`, canonical `raw_results.*`, `sse_results.csv`, and
`sse_summary.rds`. By default, replicate datasets and fitted objects are also
retained under `simulations/` and `fits/`; control those artifacts with
`saveDatasets` and `saveFits`.

Rerunning using the same output directory with `restart = FALSE` resumes or
reloads the existing analysis. Use
`runSSEControl(addModels = TRUE)` to fit only new alternatives on saved
datasets from a completed run, `recomputeSSE()` to rebuild `sse_results` /
`sse_summary` outputs from saved raw results without rerunning simulations or
fits, and `simulationPostProcess` when you need to tweak each simulated data
set before estimation. By default, recomputation writes numbered
`sse_results_recompute<k>.csv` and `sse_summary_recompute<k>.rds` files; use
`overwrite = TRUE` only when the primary summary files should be replaced.

## Plotting

`nlmixr2sse` includes plotting helpers that return `ggplot2` or combined
`patchwork` objects (in the calls below, `sse` is an output object).

**`plotSSEParameterBias(sse)`** — Bar chart of a summary statistic (default:
relative bias, %) for each parameter × model combination, faceted by parameter.
Use this to check whether a model recovers its parameters without systematic
over- or under-estimation. Relative bias near zero is what we want; large bars
indicate the presence of identifiability or model-misspecification problems.

**`plotSSEParameterEstimates(sse)`** — Boxplots and jittered points of the
replicate-level parameter estimates, faceted by parameter, with optional
overlays of the true parameter values (× marks). Use this alongside the
bias plot to distinguish a small mean bias (good) from a wide, noisy
distribution (low precision). A model can be unbiased on
average but too imprecise to be helpful.

**`plotSSEOfvDistribution(sse)`** — ECDF or histogram of
Δ OFV = OFV(reference) − OFV(alternative) across replicates. Under the null
(data generated by the simpler model), ΔOFV should cluster near zero. Positive
values mean the alternative has the lower OFV and fits better; negative values
mean the reference has the lower OFV and fits better. This helper uses the
legacy reference-versus-alternative sign convention; explicit comparisons use
the unambiguous statistic `OFV(reduced) - OFV(full)` described above.

**`plotSSEPower(sse, direction = "power")`** — A legacy empirical power curve:
the percentage of replicates in which ΔOFV is less than the negative threshold.
The companion `direction = "type1"` curve shows the percentage in which ΔOFV
exceeds the positive threshold. The default thresholds are `0`,
`3.84`, `5.99`, `7.81`, `9.49`, and `10.83`; for a one-degree-of-freedom
chi-square test, the nonzero values correspond to upper-tail p-values of
approximately 0.05, 0.014, 0.0052, 0.0021, and 0.001. At a chosen threshold,
read the Y value to obtain the estimated percentage.

**`plotSSEPpePower(sse)`** — Parametric power estimation (PPE): extrapolates
power to other study sizes by scaling an estimated non-centrality parameter.
Defaults to `method = "distribution_mle"`: one non-centrality per
`sseComparison()`, fit by maximum likelihood to the whole retained
test-statistic distribution
([Ueckert, Karlsson & Hooker 2016](https://link.springer.com/article/10.1007/s10928-016-9468-y),
the method PsN uses); the shaded ribbon is a parametric-bootstrap interval under
the fitted model (`ppeSummary()` gives the same estimate as a table, and
`plotSSEPpeDiagnostics()` checks how well the fit actually describes the
data). The previous per-threshold estimator is still available as
`method = "exceedance"`, whose ribbon is a Monte Carlo confidence interval
from the binomial uncertainty around the observed detection rate. Use either
to read off the sample size required to reach a target power (e.g., 80%). The
distribution-MLE interval is model-based: it describes estimator variability
assuming the fitted noncentral chi-square model is correct, not uncertainty
from model misspecification. Non-positive test statistics are outside that
model's support and are reported and excluded according to the selected
non-positive-statistic policy.

For PPE, use an SSE run whose explicit comparison has a scientifically
justified `df`, rather than a custom `criticalValue`. Use the tabular and
graphical diagnostics alongside the curve:

```r
ppeSummary(ppe_sse)
plotSSEPpeDiagnostics(ppe_sse)
```

PPE extrapolation additionally assumes that noncentrality scales linearly with
study size. When comparable SSE runs are available at two or more study sizes,
check that assumption explicitly:

```r
validateSSEPpeScaling(
  list(ppe_sse_small, ppe_sse_large),
  comparisons = ppe_comparison
)
```

**`plot(sse, type = "diagnostics")`** — Side-by-side four-panel view combining
parameter bias, replicate estimates, OFV distribution, and power curve. Useful
for a quick overall health check of the SSE.

```r
plotSSEParameterBias(sse)
plotSSEParameterEstimates(sse)
plotSSEOfvDistribution(sse)
plotSSEPower(sse)
plotSSEPpePower(ppe_sse)
plotSSEPpeDiagnostics(ppe_sse)
plot(sse, type = "diagnostics")
```

## Parallel execution

`runSSEControl(workers = x)` sets how many replicates run in parallel, and
`rxThreads` sets how many `rxode2` threads each worker may use. The default,
`rxThreads = "auto"`, divides the machine's cores among the workers.

```r
sse <- runSSE(
  fit,
  samples = 500,
  seed = 42,
  control = runSSEControl(workers = 4, rxThreads = "auto")
)
```

**Note:** `rxode2`'s thread count affects simulated values even under a fixed seed, so
`rxThreads` is part of a run's reproducibility profile. The resolved count is
stored in `run_info` and checked on resume, in cases in which a mismatch aborts.
Extending a run with `addModels` reuses the saved datasets unchanged, so a mismatch there
only warns and the original recorded value is kept. To reproduce a study on
different hardware, pass the recorded integer explicitly rather than relying on
`"auto"`.

## Credit where it's due

`nlmixr2sse` is based in part on the
[PsN implementation](https://github.com/UUPharmacometrics/PsN/releases/download/v5.7.0/sse_userguide.pdf).
PsN is written by Lars Lindbom,
Niclas Jonsson, Pontus Pihlgren, Mats Karlsson, Andrew Hooker, Kajsa Harling, 
Rikard Nordgren and Svetlana Freiberga.
