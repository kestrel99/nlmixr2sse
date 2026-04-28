# nlmixr2sse 0.1

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
