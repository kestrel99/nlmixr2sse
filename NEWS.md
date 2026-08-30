# nlmixr2sse 0.1

* New `validateSSEPpeScaling()`: checks the proportional-noncentrality
  assumption that PsN-style parametric power estimation relies on when
  extrapolating power to unstudied sample sizes. Fits the same power
  comparison across two or more completed SSE runs that differ only in study
  size, reports the fitted noncentrality-per-subject ratio for each, and
  warns when that ratio spreads over more than a tolerance (default 25%) of
  its mean -- evidence the noncentrality does not scale linearly with study
  size over the range tested. With three or more runs, also reports a
  through-origin lack-of-fit statistic. Restricted to power comparisons
  (Type-I comparisons have no noncentrality to scale) and refuses to pool
  runs whose resolved comparisons disagree (differing `df` or model labels)
  rather than silently averaging across them.

* New `plotSSEPpeDiagnostics()`: a distribution-adequacy diagnostic for the
  `distribution_mle` noncentral chi-square fit that `ppeSummary()`/
  `plotSSEPpePower()` rest on. Computes a Cramer-von Mises discrepancy
  between the empirical CDF of the retained test statistics and the FITTED
  chi-square CDF, with a parametric-bootstrap p-value (draw synthetic data
  from the fitted model, refit, recompute the discrepancy, and see how often
  a correctly-specified world produces one at least this large), and renders
  an ECDF-vs-fitted-CDF panel with a pointwise bootstrap envelope alongside a
  QQ panel (combined via `patchwork`, same as `plotSSEDiagnostics()`). Both
  `power` and `type1` comparisons are diagnosable -- the diagnosed parameter
  is `ncp` for a power comparison, `df` for a Type-I one, mirroring
  `.ppeParametricBootstrap()`'s existing `target` split. Multiple comparisons
  facet the plot and add one row each to the `"ppeDiagnostics"` attribute
  (`comparison`, `cvm`, `p_value`, `n`, `n_nonpositive`, `df`, `df_source`),
  so the numeric diagnostics never require reading pixels. This is evidence
  about approximation adequacy, not a pass/fail certification: a small
  `p_value` says the fitted distribution poorly describes the data, but does
  not by itself invalidate the point estimate, and a large one does not
  prove the assumption exactly right -- only not detectably wrong at this
  sample size. A comparison built with `criticalValue` instead of `df` is
  refused, same as `ppeSummary()`.
* New minimal `_pkgdown.yml` with a `reference:` index grouping the exported
  functions (run/configure, comparisons, PPE, plotting).
* New `ppeSummary()`: one row per comparison combining the `distribution_mle`
  point estimate with a parametric-bootstrap uncertainty interval and
  provenance (`df_source`, `boundary`, retained/discarded counts). The
  bootstrap draws synthetic test statistics from the FITTED noncentral
  chi-square, refits the MLE to each draw, and takes percentile quantiles of
  the refits -- always reported as `interval_type = "model_based"`, never
  "empirical": `rchisq()` draws are strictly positive and can never reproduce
  the truncation the real data underwent, so the interval covers estimator
  variability under the fitted model only, never model misspecification. This
  also fills in `plotSSEPpePower()`'s `power_lower`/`power_upper`, previously
  `NA` under `distribution_mle`, with the same bootstrap. `bootstrapSamples =
  0` skips the bootstrap (in either function) and returns the point estimate
  and counts with `NA` bounds, for callers who want the cheap estimate only.
* `plotSSEPpePower()` now renders Type-I comparisons (where the simulation
  model is the *reduced* member) as a point-range of the estimated Type-I
  rate against a dashed nominal `alpha` reference line, faceted by threshold
  like the power curve -- a Type-I comparison estimates `df`, not a
  noncentrality, so it has no sample-size curve to draw. A call that mixes
  power and Type-I comparisons renders only the power curve and warns that
  the Type-I ones were omitted, pointing to `ppeSummary()` for their
  estimates.
* `plotSSEPpePower()`/`.ppePowerPlotData()` gain `method = c("distribution_mle",
  "exceedance")`. **`"distribution_mle"` is the new default**: it fits ONE
  noncentrality parameter to the whole retained likelihood-ratio test-statistic
  distribution by maximum likelihood (Ueckert, Karlsson & Hooker 2016, the
  method PsN implements) and scales it linearly with study size, replacing
  the previous per-threshold exceedance inversion, which used only one bit
  of information per threshold and could imply a different effect size at
  every threshold. The old estimator remains available as `method =
  "exceedance"`, unchanged in behaviour, with its previously-internal
  per-threshold values now exposed as `threshold_exceedance_probability` and
  `threshold_implied_ncp` so the two estimators' outputs can never be
  confused. `distribution_mle`'s `power_lower`/`power_upper` were `NA` when
  this was first added (no ribbon was drawn); see the parametric-bootstrap
  entry above -- they are now real. PPE now requires the comparison to be built
  with `sseComparison(df = )`, not `sseComparison(criticalValue = )`:
  `distribution_mle` assumes a noncentral chi-square alternative that an
  explicit `criticalValue` gives no basis for, and refuses such comparisons
  outright.
* `runSSEControl(randomEstimationInits = )` is soft-deprecated in favor of
  `referenceInitials`/`alternativeInitials`, which split the same
  starting-value policy by model role: `referenceInitials` controls the
  simulation (reference) model refit, `alternativeInitials` controls
  alternative-model refits, so a reference-only or alternative-only
  sensitivity study is now possible. `randomEstimationInits = TRUE`/`FALSE`
  keeps working (with a deprecation warning) and maps to
  `referenceInitials`/`alternativeInitials` `"simulation"`/`"model"` for both
  roles. Setting either new argument to `"simulation"` now warns, rather than
  silently doing nothing, when `parameterSource` is not `"rawres"` -- the only
  mode the setting currently affects. A resumed run aborts if its recorded
  `referenceInitials`/`alternativeInitials` do not match the current call,
  the same way a `parameterSource` or `rxThreads` mismatch does.
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
* `runSSEControl()` gains `rxThreads`, which sets how many rxode2 threads each
  parallel worker may use. It defaults to `"auto"`, dividing the machine's
  cores among the workers. Parallel SSE runs previously aborted at worker
  counts where the workers' combined rxode2 threads exceeded the core count,
  and the setting needed to avoid that was not reachable through
  `runSSEControl()`. Requires `nlmixr2utils` >= 0.3.
* **Reproducibility note:** rxode2's thread count changes simulated values even
  under a fixed seed, so `rxThreads` affects results, not only speed. The
  resolved thread count is now recorded in `run_info` and checked when a run is
  resumed, where a mismatch aborts because the replicates would no longer be
  comparable. Extending a completed run with `addModels` refits against the
  saved simulated datasets, so a mismatch there only warns and the originally
  recorded value is kept. Because the default is `"auto"`, which derives from the machine's
  core count, reproducing a study on different hardware requires passing the
  recorded integer explicitly, e.g. `runSSEControl(rxThreads = 16)`. Runs made
  with earlier versions used rxode2's own default and will not reproduce
  byte-for-byte under the new default.
* Phase 8 adds `plotSSEPpePower()` plus persisted study-size metadata in
  `runInfo`, allowing PsN-style parametric power estimation curves with
  Monte-Carlo uncertainty ribbons to be derived from a single SSE sample size.
* Phase 7 adds `simulationPostProcess` hooks for per-replicate simulated-data
  rewriting before estimation, a plotting API built around tidy SSE outputs
  (`plotSSEParameterBias()`, `plotSSEParameterEstimates()`,
  `plotSSEOfvDistribution()`, `plotSSEPower()`, and `plotSSEDiagnostics()`),
  plus the `plot.nlmixr2SSE()` dispatcher.
* Phase 6 adds resumable `sse_state.rds` tracking, `runSSEControl(addModels = TRUE)` for fitting new alternatives onto completed runs with saved datasets, and `recomputeSSE()` for rebuilding `sse_results` / `sse_summary` outputs from saved raw results without rerunning simulations or fits.
* `runSSE()` now supports the Phase 5 uncertainty modes: canonical `rawresInput` with `inFilter`/`outFilter`, covariance-driven theta draws, `randomEstimationInits`, and `updateFix`, while continuing to write `initial_estimates.csv`, canonical `raw_results.*`, and summary outputs.
* `print()` and `summary()` methods for `nlmixr2SSE` now provide wide, human-readable summary tables on top of the tidy summary data stored in the returned object.
