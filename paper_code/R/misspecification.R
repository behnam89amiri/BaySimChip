# =============================================================================
# R/misspecification.R - virtual "truth" data and refits with an incomplete
# analysis model (Figs. 4C, 5C, 6C and 7)
#
# Model: NSBpara2DAnalysis.rds + apical solid compartment (6 compartments).
# Only the liver compartment (3) is "observed". CL_int (compartment 3) and PS
# (gut barrier, PermData[[1]]) are estimated with AnalysisODE().
# =============================================================================

IDX_APICAL <- 1
IDX_LIVER  <- 3
IDX_PS_GUT <- 1   # PermData index of the gut barrier
IDX_PS_NSB <- 2   # PermData index of the basolateral NSB surface

# Base model of the misspecification analyses (statement order as in the
# original analysis).
#   kev_apical, kev_liver : evaporation rates (mL/h)
#   ps_nsb                : basolateral PS_NSB (mL/h); apical = ps_nsb / nsb_ratio
#   k_media_plastic       : K_media/plastic of both solid compartments
#   samp_vol              : sample volume (mL)
#   times_liver/apical    : sampling schedules; the last listed time is not
#                           applied as a withdrawal (as in the original analysis) and the
#                           t = 0 withdrawal is ignored by the engine
#   final_time            : simulated (= observed) horizon (h)
#   CL, PS                : placeholder compound (overwritten by the callers)
build_misspec_app <- function(model_file,
                              kev_apical, kev_liver,
                              ps_nsb, k_media_plastic = 1, nsb_ratio = 5,
                              samp_vol = 0.05,
                              times_liver = c(0, 1, 2, 4, 8, 24, 48, 72),
                              times_apical = c(24, 48, 72),
                              final_time = 72,
                              volumes = c(0.35, 0.8, 1.6, 0.2, 1, 1),
                              CL = 0.15, PS = 0.0015,
                              molar_mass = 300, dose_uM = 1) {
  app <- load_model(model_file)
  app <- add_apical_nsb(app)
  app$nsbRatio <- nsb_ratio
  app$FinalTime <- final_time

  app$CompData[[IDX_APICAL]]$Kev <- kev_apical
  app$CompData[[IDX_LIVER]]$Kev  <- kev_liver

  app$CompData[[5]]$fu <- k_media_plastic
  app$CompData[[6]]$fu <- k_media_plastic
  app$PermData[[2]]$PS <- ps_nsb
  app$PermData[[3]]$PS <- app$PermData[[2]]$PS / app$nsbRatio

  for (i in seq_along(volumes)) app$CompData[[i]]$Vol <- volumes[i]

  app$CompData[[IDX_LIVER]]$ExMedia  <- exchange_table(times_liver, samp_vol)
  app$CompData[[IDX_APICAL]]$ExMedia <- exchange_table(times_apical, samp_vol)

  app$PermData[[IDX_PS_GUT]]$PS <- PS
  app$CompData[[IDX_LIVER]]$CL  <- CL

  app$MolarMass <- molar_mass
  app$CompData[[1]]$DoseCond <- 1
  app$CompData[[1]]$DoseSingData <- app$CompData[[1]]$Vol * dose_uM * app$MolarMass

  set_fixed_horizon_and_expdata(app, Tfin = app$FinalTime, times_liver = times_liver)
}

# Observation times (liver only) and horizon; trims media-exchange tables to
# events <= Tfin
set_fixed_horizon_and_expdata <- function(app, Tfin, times_liver) {
  app$FinalTime <- Tfin
  if (is.null(app$expData)) app$expData <- list()
  while (length(app$expData) < 3) app$expData[[length(app$expData) + 1]] <- list()
  app$expData[[1]] <- list(t = NULL, y = NULL)
  app$expData[[2]] <- list(t = NULL, y = NULL)
  app$expData[[3]] <- list(t = sort(unique(c(times_liver[times_liver <= Tfin], Tfin))), y = NULL)
  for (i in seq_along(app$CompData)) {
    ex <- app$CompData[[i]]$ExMedia
    if (!(is.null(ex) || !is.matrix(ex) || ncol(ex) < 3 || nrow(ex) == 0)) {
      app$CompData[[i]]$ExMedia <- ex[ex[, 1] <= Tfin, , drop = FALSE]
    }
  }
  app
}

# Median non-zero sample volume of a compartment (NA if none)
get_base_sampling_volume <- function(app, comp_idx) {
  ex <- app$CompData[[comp_idx]]$ExMedia
  if (is.null(ex) || ncol(ex) < 2 || nrow(ex) == 0) return(NA_real_)
  v <- suppressWarnings(as.numeric(ex[, 2]))
  v <- v[is.finite(v) & v > 0]
  if (length(v) == 0) return(NA_real_)
  stats::median(v)
}

# Replace the sample volume of a compartment; volume 0 removes all events
set_sampling_volume <- function(app, comp_idx, v_mL) {
  empty_exmedia <- matrix(numeric(0), ncol = 3)
  if (!is.finite(v_mL) || v_mL <= 0) {
    app$CompData[[comp_idx]]$ExMedia <- empty_exmedia
    return(app)
  }
  ex <- app$CompData[[comp_idx]]$ExMedia
  if (!is.null(ex) && nrow(ex) > 0 && ncol(ex) >= 1) {
    times <- as.numeric(ex[, 1])
  } else {
    times <- tryCatch(as.numeric(app$expData[[comp_idx]]$t), error = function(e) numeric(0))
  }
  if (length(times) == 0) {
    app$CompData[[comp_idx]]$ExMedia <- empty_exmedia
  } else {
    app$CompData[[comp_idx]]$ExMedia <- cbind(times, rep(v_mL, length(times)), rep(0, length(times)))
  }
  app
}

# Simulate the reference ("truth") model and sample the liver compartment at
# the observation times, with proportional noise: y = max(C x (1 + e), 0),
# e ~ N(0, prop_noise_sd). Draws one rnorm() per observation, in time order.
simulate_truth <- function(app_in, prop_noise_sd = 0) {
  app <- copy_app(app_in)
  app <- ConstructDoseMediaEx(app)
  app <- ConstructODE(app)
  app <- CalCulate(app)
  times_needed <- app$expData[[IDX_LIVER]]$t
  app$expData[[IDX_LIVER]]$y <- sapply(times_needed, function(tt) {
    i <- which.min(abs(app$t - tt))
    val <- app$yc[i, IDX_LIVER]
    if (prop_noise_sd > 0) val <- max(val * (1 + stats::rnorm(1, 0, prop_noise_sd)), 0)
    val
  })
  app$expData[[1]]$y <- NULL
  app$expData[[2]]$y <- NULL
  app
}

pct_bias <- function(est, truth) 100 * (est - truth) / truth

# bias (%) -> absolute fold error (+100% -> 2x, -50% -> 2x); NA if bias <= -100%
bias_to_fold <- function(bias_pct) {
  x <- 1 + bias_pct / 100
  out <- rep(NA_real_, length(x))
  ok <- is.finite(x) & x > 0
  out[ok] <- pmax(x[ok], 1 / x[ok])
  out
}

# Estimate CL_int and PS from the truth data with an analysis model in which
# one feature is set to `value` (Fig. 7A):
#   "evaporation": value = list(apical = Kev, liver = Kev)
#   "sampling"   : value = list(apical = Vsamp, liver = Vsamp)
#   "nsb"        : value = PS_NSB (basolateral; apical = value / nsbRatio)
fit_with_misspecification <- function(app_truth, feature, value, r0) {
  app_fit <- copy_app(app_truth)
  if (feature == "evaporation") {
    app_fit$CompData[[IDX_APICAL]]$Kev <- value$apical
    app_fit$CompData[[IDX_LIVER]]$Kev  <- value$liver
  } else if (feature == "sampling") {
    app_fit <- set_sampling_volume(app_fit, IDX_APICAL, value$apical)
    app_fit <- set_sampling_volume(app_fit, IDX_LIVER, value$liver)
  } else if (feature == "nsb") {
    app_fit$PermData[[IDX_PS_NSB]]$PS <- value
    if (length(app_fit$PermData) >= 3 && is.finite(app_fit$nsbRatio) && app_fit$nsbRatio > 0) {
      app_fit$PermData[[3]]$PS <- value / app_fit$nsbRatio
    }
  }
  app_fit$Estimated <- list(list(Done = 0), list(Done = 0))
  app_fit$r0 <- r0
  app_fit$CompData[[IDX_LIVER]]$CL <- NULL
  app_fit$PermData[[IDX_PS_GUT]]$PS <- NULL
  app_fit <- ConstructDoseMediaEx(app_fit)
  app_fit <- ConstructODE(app_fit)
  app_fit <- AnalysisODE(app_fit)
  list(est_CL = app_fit$Estimated[[1]]$Values, est_PS1 = app_fit$Estimated[[2]]$Values)
}

# One compound: truth with all processes at nominal values, then refits with
# each process at level x nominal (level 0 = process neglected). Returns one
# row per process and level for this compound.
run_all_sweeps <- function(app_base, levels, prop_noise_sd) {
  base_kev_apical   <- app_base$CompData[[IDX_APICAL]]$Kev
  base_kev_liver    <- app_base$CompData[[IDX_LIVER]]$Kev
  base_vsamp_apical <- get_base_sampling_volume(app_base, IDX_APICAL)
  base_vsamp_liver  <- get_base_sampling_volume(app_base, IDX_LIVER)
  base_ps_nsb       <- app_base$PermData[[IDX_PS_NSB]]$PS
  bases <- c(base_kev_apical, base_kev_liver, base_vsamp_apical, base_vsamp_liver, base_ps_nsb)
  if (any(!is.finite(bases) | bases <= 0)) stop("All nominal process parameters must be > 0.")

  kev_apical_grid   <- levels * base_kev_apical
  kev_liver_grid    <- levels * base_kev_liver
  vsamp_apical_grid <- levels * base_vsamp_apical
  vsamp_liver_grid  <- levels * base_vsamp_liver
  ps_nsb_grid       <- levels * base_ps_nsb

  truth_CL  <- app_base$CompData[[IDX_LIVER]]$CL
  truth_PS1 <- app_base$PermData[[IDX_PS_GUT]]$PS
  app_truth <- simulate_truth(app_base, prop_noise_sd = prop_noise_sd)
  r0 <- c(app_truth$CompData[[IDX_LIVER]]$CL, app_truth$PermData[[IDX_PS_GUT]]$PS)

  results <- list()
  for (k in seq_along(levels)) {
    fit <- fit_with_misspecification(app_truth, "evaporation",
                                     list(apical = kev_apical_grid[k], liver = kev_liver_grid[k]), r0)
    results[[length(results) + 1]] <- data.frame(
      scenario = "mispec_evaporation", varied_param = "Kev_apical+Kev_liver", value = levels[k],
      val_apical = kev_apical_grid[k], val_liver = kev_liver_grid[k],
      est_CL = fit$est_CL, est_PS1 = fit$est_PS1,
      bias_CL_pct = pct_bias(fit$est_CL, truth_CL), bias_PS1_pct = pct_bias(fit$est_PS1, truth_PS1))
  }
  for (k in seq_along(levels)) {
    fit <- fit_with_misspecification(app_truth, "sampling",
                                     list(apical = vsamp_apical_grid[k], liver = vsamp_liver_grid[k]), r0)
    results[[length(results) + 1]] <- data.frame(
      scenario = "mispec_sampling", varied_param = "Vsamp_apical+Vsamp_liver", value = levels[k],
      val_apical = vsamp_apical_grid[k], val_liver = vsamp_liver_grid[k],
      est_CL = fit$est_CL, est_PS1 = fit$est_PS1,
      bias_CL_pct = pct_bias(fit$est_CL, truth_CL), bias_PS1_pct = pct_bias(fit$est_PS1, truth_PS1))
  }
  for (k in seq_along(levels)) {
    fit <- fit_with_misspecification(app_truth, "nsb", ps_nsb_grid[k], r0)
    results[[length(results) + 1]] <- data.frame(
      scenario = "mispec_nsb", varied_param = "PS_nsb", value = levels[k],
      val_nsb = ps_nsb_grid[k],
      est_CL = fit$est_CL, est_PS1 = fit$est_PS1,
      bias_CL_pct = pct_bias(fit$est_CL, truth_CL), bias_PS1_pct = pct_bias(fit$est_PS1, truth_PS1))
  }
  dplyr::bind_rows(results)
}

# All archetypes of Fig. 7A: archetypes = named list
# of c(CL_int, PS) in mL/h; the RNG stream continues across archetypes.
run_archetypes <- function(app, archetypes, levels, prop_noise_sd) {
  parts <- list()
  for (nm in names(archetypes)) {
    app$PermData[[IDX_PS_GUT]]$PS <- archetypes[[nm]][2]
    app$CompData[[IDX_LIVER]]$CL  <- archetypes[[nm]][1]
    res <- run_all_sweeps(app, levels, prop_noise_sd)
    res$name <- nm
    parts[[nm]] <- res
  }
  dplyr::bind_rows(parts)
}

# Fold errors of the neglected processes (level 0) per archetype and process
fold_errors <- function(results) {
  axis_map <- c(mispec_evaporation = "Evaporation", mispec_sampling = "Sampling", mispec_nsb = "NSB")
  r <- results[results$value == 0, ]
  data.frame(name = r$name, axis = unname(axis_map[r$scenario]),
             fold_CL = bias_to_fold(r$bias_CL_pct), fold_PS = bias_to_fold(r$bias_PS1_pct))
}

# Alternative Fig. 7B analysis (run_alldrugs(analysis = "apical_nsb_kept")): refit
# of one drug with evaporation, sampling and the basolateral NSB surface removed
# from the analysis model; the apical NSB surface is kept (PS_NSB / nsbRatio of
# the truth model, K of the apical solid compartment unchanged).
fit_alldrugs_neglected <- function(app_truth, app) {
  app_fit <- copy_app(app_truth)
  app_fit$CompData[[1]]$Kev <- 0
  app_fit$CompData[[3]]$Kev <- 0
  app_fit$CompData[[5]]$fu <- 1
  app_fit$CompData[[3]]$fu <- 1
  app_fit$PermData[[2]]$PS <- 0
  app_fit$PermData[[3]]$PS <- app$PermData[[2]]$PS / app$nsbRatio
  app_fit$CompData[[1]]$ExMedia[, 2] <- 0
  app_fit$CompData[[3]]$ExMedia[, 2] <- 0
  app_fit$Estimated <- list(list(Done = 0), list(Done = 0))
  app_fit$r0 <- c(app_truth$CompData[[3]]$CL, app_truth$PermData[[1]]$PS)
  app_fit$CompData[[3]]$CL <- NULL
  app_fit$PermData[[1]]$PS <- NULL
  app_fit <- ConstructDoseMediaEx(app_fit)
  app_fit <- ConstructODE(app_fit)
  app_fit <- AnalysisODE(app_fit)
  c(app_fit$Estimated[[1]]$Values, app_fit$Estimated[[2]]$Values)
}

# All drugs of Fig. 7B; drugs: data frame with Compound, CL, PS (mL/h). A
# failing fit gives NA.
#   analysis = "all_neglected"  : evaporation, sampling and NSB (both solid
#              compartments) removed from the analysis model, as for the bias
#              maps (fit_neglecting_all()); used for Fig. 7B
#   analysis = "apical_nsb_kept": evaporation, sampling and basolateral NSB
#              removed, apical NSB kept (fit_alldrugs_neglected())
run_alldrugs <- function(app, drugs, prop_noise_sd, analysis = c("all_neglected", "apical_nsb_kept")) {
  analysis <- match.arg(analysis)
  results <- list()
  for (i in seq_len(nrow(drugs))) {
    app$PermData[[1]]$PS <- drugs$PS[i]
    app$CompData[[3]]$CL <- drugs$CL[i]
    truth_CL  <- app$CompData[[3]]$CL
    truth_PS1 <- app$PermData[[1]]$PS
    app_truth <- simulate_truth(app, prop_noise_sd = prop_noise_sd)
    est <- tryCatch(if (analysis == "apical_nsb_kept") fit_alldrugs_neglected(app_truth, app)
                    else unname(fit_neglecting_all(app_truth, truth_CL, truth_PS1)), error = function(e) {
      warning(drugs$Compound[i], ": fit failed (", conditionMessage(e), ")")
      c(NA_real_, NA_real_)
    })
    results[[i]] <- data.frame(
      Compound = drugs$Compound[i], scenario = "All_phenomena_neglected",
      est_CL_raw = est[1], est_PS1_raw = est[2],
      bias_CL_pct = pct_bias(est[1], truth_CL), bias_PS1_pct = pct_bias(est[2], truth_PS1),
      est_CL = est[1] * 10 / 6 * 1000 / 60,                 # µL/min/10^6 hepatocytes
      est_P = est[2] * 1000000 / (3600 * GUT_AREA_CM2),       # 10^-6 cm/s
      true_CL = truth_CL * 10 / 6 * 1000 / 60,
      true_P = truth_PS1 * 1000000 / (3600 * GUT_AREA_CM2))
  }
  dplyr::bind_rows(results)
}

# Estimate CL_int and PS with evaporation, sampling and NSB all removed from
# the analysis model (bias maps and Fig. 7B)
fit_neglecting_all <- function(app_truth, CL_true, PS_true) {
  app_fit <- copy_app(app_truth)
  app_fit$CompData[[IDX_APICAL]]$Kev <- 0
  app_fit$CompData[[IDX_LIVER]]$Kev  <- 0
  for (idx in c(IDX_APICAL, IDX_LIVER)) {
    ex <- app_fit$CompData[[idx]]$ExMedia
    if (!is.null(ex) && nrow(ex) > 0) app_fit$CompData[[idx]]$ExMedia[, 2] <- 0
  }
  if (length(app_fit$CompData) >= 5) app_fit$CompData[[5]]$fu <- 1
  if (length(app_fit$CompData) >= 6) app_fit$CompData[[6]]$fu <- 1
  app_fit$CompData[[IDX_LIVER]]$fu <- 1
  if (length(app_fit$PermData) >= IDX_PS_NSB) app_fit$PermData[[IDX_PS_NSB]]$PS <- 0
  if (length(app_fit$PermData) >= 3) app_fit$PermData[[3]]$PS <- 0

  app_fit$Estimated <- list(list(Done = 0), list(Done = 0))
  app_fit$r0 <- c(CL_true, PS_true)
  app_fit$CompData[[IDX_LIVER]]$CL <- NULL
  app_fit$PermData[[IDX_PS_GUT]]$PS <- NULL
  app_fit <- ConstructDoseMediaEx(app_fit)
  app_fit <- ConstructODE(app_fit)
  app_fit <- AnalysisODE(app_fit)
  c(est_CL = app_fit$Estimated[[1]]$Values, est_PS = app_fit$Estimated[[2]]$Values)
}
