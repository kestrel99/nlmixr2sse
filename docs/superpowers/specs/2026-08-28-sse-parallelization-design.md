# SSE Parallelization Design

## Goal

Support the `nlmixr2utils` 0.3 parallelization API in `nlmixr2sse`: expose
`rxThreads`, thread it through the worker-plan and task-apply calls, and record
it as the reproducibility-relevant setting it turns out to be.

Parallel SSE is currently unusable at any meaningful worker count. This is a
defect fix, independent of the two feature specs dated the same day.

## Background

### The defect

`nlmixr2utils` 0.3 added a thread-oversubscription guard to `.withWorkerPlan()`
and documented it as a breaking change:

> existing calls to `.withWorkerPlan()` with `workers > 1` or `workers = "auto"`
> that do not also set `rxThreads` will, on most real multi-core machines, now
> abort where they previously ran silently oversubscribed.

`nlmixr2sse` makes exactly those calls. Both call sites
(`run-sse.R:412`, `run-sse.R:480`) pass no `rxThreads`, and neither do the
inner `.plap()` calls. On a 32-core machine with rxode2 defaulting to 16
threads:

```
workers = 2  ->  32 threads vs 32 cores  ->  runs (at the limit)
workers = 4  ->  64 threads vs 32 cores  ->  ABORTS
```

```
Requested 4 workers x 16 rxode2 threads = 64 threads, but only 32 cores
are available.
Lower `workers`, or set `rxThreads` (e.g. `rxThreads = 1`) so their product is <= 32.
```

`runSSEControl()` has **no `rxThreads` argument**, so the remedy the error
message names is unreachable through the package's public API. The only
workaround is to lower `workers` — on this machine, to 2.

For a tool whose normal usage is hundreds to a thousand replicates, losing
parallel execution is a serious limitation.

### Thread count changes results

Thread count is not merely a performance knob. rxode2 distributes subjects
across threads with per-thread RNG streams, so the number of threads changes
the simulated data even with a fixed seed.

Verified directly — same model, same events, same `set.seed()` and
`rxSetSeed()`, varying only `setRxThreads()`:

```
nthreads 1 vs 4 identical: FALSE   max abs diff: 0.1817
nthreads 1 vs 8 identical: FALSE   max abs diff: 0.1491
```

`rxThreads` is therefore part of a run's reproducibility contract and must be
recorded, validated on resume, and documented. This hazard already exists
today — rxode2's own default derives from the machine's core count — but it is
currently invisible and unrecorded.

## Chosen Approach

Add `rxThreads` to `runSSEControl()`, pass it to both `.withWorkerPlan()` calls
and both `.plap()` calls, and persist the **resolved** value in `run_info`.

### Default

`rxThreads = "auto"`.

`nlmixr2utils::resolveRxThreads()` divides available cores among workers:

| `workers` | `"auto"` resolves to |
| --- | --- |
| 1 | 32 |
| 2 | 16 |
| 4 | 8 |
| 8 | 4 |
| `"auto"` | 1 |
| `NULL` | 32 |

(on a 32-core machine)

This fixes the abort at every worker count and uses the machine better than
today's fixed 16 threads for single-worker runs.

**The trade-off, accepted deliberately:** because `"auto"` derives from core
count, the same call with the same seed produces different simulated data on
machines with different core counts. This is mitigated, not eliminated:

- The **resolved** integer — not the string `"auto"` — is written to
  `run_info`, so any completed run states exactly what it used.
- The run banner already prints the resolved thread count per worker.
- Documentation states that reproducing a study on different hardware requires
  setting `rxThreads` to the recorded integer explicitly.

Note that today's `NULL` default is *also* machine-dependent, so `"auto"` does
not introduce the hazard — it makes it visible and recoverable.

## Component Design

| Function | Change |
| --- | --- |
| `runSSEControl()` | New `rxThreads = "auto"` argument, validated as a positive integer or `"auto"` (delegating to `nlmixr2utils`' own validator where possible). Stored on the control object. |
| `runSSE()` | Resolve `rxThreads` once via `nlmixr2utils::resolveRxThreads(control$workers, control$rxThreads)`; record the resolved integer in `run_info`; pass `rxThreads` to both `.withWorkerPlan()` calls and both `.plap()` calls. |
| `.workerDescription()` | Report the resolved thread count alongside the worker description, so the SSE banner states the full parallel configuration in one line. |
| `.validateResumeRequest()` | Compare the recorded resolved `rxThreads` against the current run's resolved value; abort on mismatch (see below). |
| `DESCRIPTION` | Raise `nlmixr2utils` from `>= 0.2` to `>= 0.3` — `resolveRxThreads()` and the `rxThreads` arguments are 0.3 features. |

Resolution happens **once** in `runSSE()`, and the resolved integer is what is
passed onward and recorded. Resolving separately at each call site could yield
different values if the ambient `future` plan changed mid-run.

### Resume behaviour

Because thread count changes simulated data, resuming a partially complete run
under a different `rxThreads` would produce replicates drawn under two
different configurations within one study.

On resume, compare the recorded resolved value against the current resolved
value and abort on mismatch, naming both and stating that the original value
can be restored with `runSSEControl(rxThreads = <recorded>)`. This mirrors the
existing resume validation for `samples`, `parameterSource`,
`randomEstimationInits`, and `updateFix`.

Runs resumed from directories written before this change carry no recorded
value. Treat a missing record as "unknown" and warn rather than abort — the
original run's thread count is genuinely unknowable, so refusing to resume
would strand existing work while a warning states the risk plainly.

## Error Handling

| Condition | Behaviour |
| --- | --- |
| `rxThreads` not a positive integer or `"auto"` | Abort from `runSSEControl()` at construction. |
| `workers` x `rxThreads` exceeds core count | `nlmixr2utils` aborts with its own message naming both; not re-wrapped. |
| Resume with a different resolved `rxThreads` | Abort naming recorded and current values, and how to restore. |
| Resume from a directory with no recorded value | Warn that the original thread count is unknown and results may not be comparable; proceed. |

## Testing

Unit tests:

- `runSSEControl()` accepts `"auto"`, accepts a positive integer, and rejects
  `0`, negatives, non-integers, and arbitrary strings.
- The default is `"auto"`.
- The resolved integer, not the literal `"auto"`, is what reaches `run_info`.
- Resume with a mismatched recorded value aborts naming both values.
- Resume with no recorded value warns and proceeds.

Integration tests:

- A run with `workers = 1` completes and records a resolved thread count.
- A worker count that previously aborted (`workers = 4` on a machine where
  `workers x` rxode2's default exceeds core count) now completes. Guarded by a
  core-count skip so it does not fail on small CI machines.
- `rxThreads` composes with each `parameterSource` mode.

Deliberately **not** tested: that two different thread counts produce identical
results. They do not, by design of rxode2's threading, and asserting otherwise
would encode a false expectation.

## Risks

- **Cross-machine reproducibility.** Accepted and mitigated as described. The
  residual risk is a user who reruns a study on different hardware and does not
  consult `run_info`. Documentation and the run banner are the mitigations.
- **Behaviour change for single-worker runs.** The default moves from rxode2's
  own thread count to `"auto"` (32 vs 16 on this machine), so existing
  single-worker runs will produce different simulated data than before at the
  same seed. This must be a NEWS entry; it is the same class of change as the
  parameter-uncertainty spec's reproducibility break, and the two will ship
  visible together.
- **`nlmixr2utils` dependency bump** to `>= 0.3` affects anyone pinned to 0.2.
  Unavoidable — the API being adopted does not exist in 0.2.

## Out of Scope

- The `nlmixr2utils:::.abortRawResults()` `cli` `.envir` defect, which masks
  the covariance/rawres integration test. Filed separately; it is a bug in the
  sibling package, not a parallelization concern.
- Changing rxode2's threading or RNG-stream behaviour.
- Parallelizing anything not already parallel; this spec changes how existing
  parallel sections are configured, not what runs in parallel.
