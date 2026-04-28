# Generate the precomputed artifact used by vignette("sse-power").
#
# Run line-by-line in RStudio (working directory = package root), or via:
#   Rscript inst/scripts/generate-sse-power-example.R

# ── 1. Load packages ──────────────────────────────────────────────────────────

if (!requireNamespace("nlmixr2sse", quietly = TRUE)) {
  pkgload::load_all(".", quiet = TRUE, export_all = FALSE, helpers = FALSE)
}
suppressPackageStartupMessages({
  library(nlmixr2)
  library(nlmixr2sse)
})

# ── 2. Models ─────────────────────────────────────────────────────────────────

full_model <- function() {
  ini({
    tka <- log(1.2)
    tcl <- log(4.5)
    tv <- log(45)
    tcl_wt <- 0.30
    eta.cl ~ 0.09
    eta.v ~ 0.09
    prop.sd <- 0.20
  })
  model({
    ka <- exp(tka)
    cl <- exp(tcl + tcl_wt * WTGRP + eta.cl)
    v <- exp(tv + eta.v)
    cp <- linCmt()
    cp ~ prop(prop.sd)
  })
}

base_model <- function() {
  ini({
    tka <- log(1.2)
    tcl <- log(4.5)
    tv <- log(45)
    eta.cl ~ 0.09
    eta.v ~ 0.09
    prop.sd <- 0.20
  })
  model({
    ka <- exp(tka)
    cl <- exp(tcl + eta.cl)
    v <- exp(tv + eta.v)
    cp <- linCmt()
    cp ~ prop(prop.sd)
  })
}

# ── 3. Study design ───────────────────────────────────────────────────────────

study_sizes   <- c(20L, 40L, 60L, 80L)
sample_times  <- c(0.5, 1, 2, 4, 8)
dose          <- 100
n_samples     <- 12L
seed          <- 101L

ids      <- seq_len(max(study_sizes))
wt_group <- rep(c(0L, 1L), length.out = max(study_sizes))

dose_rows <- data.frame(
  ID    = ids,
  TIME  = 0,
  DV    = NA_real_,
  AMT   = dose,
  EVID  = 101L,
  CMT   = 1L,
  WTGRP = wt_group,
  stringsAsFactors = FALSE
)
obs_rows <- expand.grid(
  ID   = ids,
  TIME = sample_times,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
obs_rows$WTGRP <- wt_group[match(obs_rows$ID, ids)]
obs_rows$DV    <- NA_real_
obs_rows$AMT   <- 0
obs_rows$EVID  <- 0L
obs_rows$CMT   <- 2L
obs_rows <- obs_rows[, names(dose_rows)]

design <- rbind(dose_rows, obs_rows)
design <- design[order(design$ID, design$TIME, -design$EVID), , drop = FALSE]
rownames(design) <- NULL

# ── 4. Simulate reference data ────────────────────────────────────────────────

set.seed(seed)

sim_ui    <- nlmixr2est::nlmixr2(full_model)
sim_theta <- c(tka = log(1.2), tcl = log(4.5), tv = log(45),
               tcl_wt = 0.30, prop.sd = 0.20)
sim_omega <- diag(c(0.09, 0.09))
dimnames(sim_omega) <- list(c("eta.cl", "eta.v"), c("eta.cl", "eta.v"))

design$.row_id <- seq_len(nrow(design))
solved <- as.data.frame(
  rxode2::rxSolve(
    sim_ui,
    params       = sim_theta,
    events       = design,
    omega        = sim_omega,
    sigma        = matrix(numeric(0), 0L, 0L),
    addCov       = TRUE,
    keep         = ".row_id",
    returnType   = "data.frame"
  ),
  stringsAsFactors = FALSE
)
solved <- solved[order(solved$.row_id), , drop = FALSE]

obs_mask             <- design$EVID == 0L
solved_obs_mask      <- solved$.row_id %in% design$.row_id[obs_mask]
design$DV[obs_mask]  <- solved$sim[solved_obs_mask]
design$.row_id       <- NULL
reference_data  <- design

# ── 5. SSE runs ───────────────────────────────────────────────────────────────

runs <- vector("list", length(study_sizes))
names(runs) <- as.character(study_sizes)

for (i in seq_along(study_sizes)) {
  n <- study_sizes[[i]]
  message("Running SSE for N = ", n)

  keep_ids   <- sort(unique(reference_data$ID))[seq_len(n)]
  study_data <- reference_data[reference_data$ID %in% keep_ids, , drop = FALSE]
  study_data$ID <- match(study_data$ID, keep_ids)
  rownames(study_data) <- NULL

  full_fit <- nlmixr2(
    full_model, study_data, est = "focei",
    control = list(print = 0L, covMethod = "r", eval.max = 80L, iter.max = 80L)
  )
  base_fit <- nlmixr2(
    base_model, study_data, est = "focei",
    control = list(print = 0L, covMethod = "r", eval.max = 80L, iter.max = 80L)
  )

  runs[[i]] <- runSSE(
    full_fit,
    alternativeModels = sseModel(base_fit, label = "no_wt"),
    samples    = n_samples,
    seed       = seed + n,
    control    = runSSEControl(workers = 1L),
    outputDir  = file.path(tempdir(), paste0("sse-power-", n)),
    restart    = TRUE
  )
}

# ── 6. Summarise results ──────────────────────────────────────────────────────

power_curve <- do.call(rbind, lapply(seq_along(runs), function(i) {
  sse        <- runs[[i]]
  study_size <- study_sizes[[i]]
  ps         <- sse$powerSummary

  power_row <- ps[
    ps$model_label == "no_wt" &
      !is.na(ps$direction) & ps$direction == "power" &
      abs(ps$threshold - 3.84) < 1e-8,
    ,
    drop = FALSE
  ]

  raw <- sse$rawResults[sse$rawResults$sample > 0L, , drop = FALSE]
  samples <- sort(unique(raw$sample))
  n_converged <- sum(vapply(samples, function(s) {
    rows     <- raw[raw$sample == s, , drop = FALSE]
    full_row <- rows[rows$role == "simulation", , drop = FALSE]
    alt_row  <- rows[rows$model_label == "no_wt", , drop = FALSE]
    nrow(full_row) == 1L && nrow(alt_row) == 1L &&
      is.na(full_row$error_message) && is.na(alt_row$error_message) &&
      is.finite(full_row$objf) && is.finite(alt_row$objf)
  }, logical(1)))

  data.frame(
    study_size  = as.integer(study_size),
    power       = as.numeric(power_row$value[[1L]]) / 100,
    n_converged = n_converged,
    stringsAsFactors = FALSE
  )
}))

effect_n   <- study_sizes[[ceiling(length(study_sizes) / 2L)]]
effect_sse <- runs[[match(effect_n, study_sizes)]]
effect_raw <- effect_sse$rawResults[effect_sse$rawResults$sample > 0L, , drop = FALSE]
samples    <- sort(unique(effect_raw$sample))

effect_estimates <- effect_raw[
  effect_raw$role == "simulation",
  c("sample", "tcl_wt"),
  drop = FALSE
]
names(effect_estimates)[2] <- "estimate"

effect_significant <- vapply(samples, function(s) {
  rows      <- effect_raw[effect_raw$sample == s, , drop = FALSE]
  full_objf <- rows$objf[rows$role == "simulation"][[1L]]
  alt_objf  <- rows$objf[rows$model_label == "no_wt"][[1L]]
  is.finite(full_objf) && is.finite(alt_objf) && (alt_objf - full_objf) > 3.84
}, logical(1))

effect_plot_data <- data.frame(
  study_size  = as.integer(effect_n),
  estimate    = as.numeric(effect_estimates$estimate),
  significant = effect_significant[match(effect_estimates$sample, samples)],
  stringsAsFactors = FALSE
)

result <- list(
  study_sizes      = as.integer(study_sizes),
  true_effect      = 0.30,
  n_for_effect_plot = as.integer(effect_n),
  power_curve      = power_curve,
  effect_plot_data = effect_plot_data
)

# ── 7. Save artifact ──────────────────────────────────────────────────────────

output_path <- file.path("inst", "extdata", "sse_power_example.rds")
saveRDS(result, output_path)
message("Saved: ", output_path)
print(power_curve)
