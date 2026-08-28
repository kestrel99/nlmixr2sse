# nlmixr2sse 0.1

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
