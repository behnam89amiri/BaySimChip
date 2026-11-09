# =============================================================================
# fig6d_empty_chip.R - Fig. 6D: NSB in an empty HUMIMIC Chip2 (PDMS),
# diclofenac 100 µM, two dosing phases followed by a washout
#
# Empty-chip NSB model: two medium
# compartments of 0.4 mL connected by the flow Q, and one solid compartment
# (reference volume 1 mL) exchanging with medium compartment 2:
#   dc1/dt = Q (c2 - c1) / 0.4
#   dc2/dt = (Q (c1 - c2) + PS_NSB (K c3 - c2)) / 0.4
#   dc3/dt = -PS_NSB (K c3 - c2) / 1          (K = K_media/plastic, "fu" in the code)
# Each dosing phase starts at the time point with the highest mean
# concentration among the first time points (incomplete mixing at early
# times); the solid compartment starts at k3 x c0 per dosing phase; the
# washout starts with medium 0 and the solid compartment at the end of the
# last dosing phase. Observed = medium compartment 1 (mean of replicates).
# Fitted: K_media/plastic, PS_NSB (and k3 if FIT_K3), least squares
# (minpack.lm::nls.lm), all conditions of both compounds.
#
# Two fit configurations (FITS) are run:
#   "joint_k3fit"       : dosing phases AND washout fitted, k3 fitted
#                         (alternative configuration; PDMS DCF 100 µM:
#                         K = 6.27, PS_NSB = 0.0234 mL/h).
#   "validate_k3fixed"  : dosing phases fitted, washout predicted, k3 = 0
#                         (solid compartment empty at the start of each dosing
#                         phase); shown in Fig. 6D. PDMS DCF 100 µM:
#                         K = 5.83, PS_NSB = 0.054 mL/h.
# PANEL_FIT selects the configuration shown in the panel.
#
# Inputs : data/empty_chip_nsb/empty_chip_nsb_dcf.xlsx, empty_chip_nsb_fiau.xlsx
# Outputs: per configuration: fig6d_fit_results_<fit>.csv, fig6d_predictions_<fit>.csv,
#          fig6d_observed_mean_<fit>.csv, fig6d_observed_raw_<fit>.csv,
#          fig6d_fit_failures_<fit>.csv; panel: fig6d_panel_*.csv, fig6d_panel.png
# Runtime: ~30 s
# =============================================================================

.f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
.root <- Sys.getenv("PAPER_CODE_ROOT", if (length(.f)) dirname(dirname(.f)) else ".")
source(file.path(.root, "R", "setup.R")); pc_setup(.root)
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(ggplot2)
})

# ---- parameters ----------------------------------------------------------------
XLSX_FILES <- c(pc_data("empty_chip_nsb", "empty_chip_nsb_dcf.xlsx"),
                pc_data("empty_chip_nsb", "empty_chip_nsb_fiau.xlsx"))
FITS <- list(
  joint_k3fit      = list(mode = "joint", fit_k3 = TRUE),
  validate_k3fixed = list(mode = "validate_washout", fit_k3 = FALSE)
)
PANEL_FIT <- "validate_k3fixed"
PANEL_CONDITION <- "PDMS__DCF__100"          # material__compound__concentration (µM)
Q_ML_H <- 2.27 * 60 / 1000                   # flow between the two medium compartments (mL/h)
V_MEDIUM_ML <- 0.4                           # each medium compartment
V_SOLID_ML  <- 1                             # reference volume of the solid compartment
PAR_INIT  <- c(fu = 0.5, PS = 0.05)          # K_media/plastic, PS_NSB (mL/h)
PAR_LOWER <- c(fu = 0.2, PS = 0.0001)
PAR_UPPER <- c(fu = 10.0, PS = 1)
K3_INIT <- 0.5; K3_LOWER <- 0; K3_UPPER <- 20   # c3(0) = k3 x c0 per dosing phase
WASHOUT_C12_INIT <- "zero"                   # medium at the start of the washout
DROP_CANDIDATES <- 1:5                       # candidate numbers of dropped early time points
DROP_WASHOUT_N <- 0
MAXITER <- 200
# excluded sheets (broken chip, mislabelled or ambiguous sheets)
EXCLUDE_SHEETS <- c("PDMS_DCF_Dose2_Conc_0uM", "nonPDMS_DCF_Dose2_Conc_     _1",
                    "PDMS_DCF_Dose2_Conc_         10", "PDMS_DCF_Dose1_Conc_        0p1",
                    "nonPDMS_FIAU_Dose1_Conc_     _1", "nonPDMS_FIAU_Dose2_Conc_     _1",
                    "PDMS_FIAU_Dose1_Conc_        10", "PDMS_FIAU_Dose2_Conc_        10")
SKIP_CONC_UM <- 0                            # controls are not fitted

# ---- model ---------------------------------------------------------------------
ode_model <- function(t, y, p) {
  c1 <- y[1]; c2 <- y[2]; c3 <- y[3]
  Q <- p["Q"]; PS <- p["PS"]; fu <- p["fu"]
  list(c(Q * (c2 - c1) / V_MEDIUM_ML,
         (Q * (c1 - c2) + PS * (fu * c3 - c2)) / V_MEDIUM_ML,
         -PS * (fu * c3 - c2) / V_SOLID_ML))
}
simulate_phase <- function(times, y0, fu, PS, Q = Q_ML_H) {
  out <- deSolve::ode(y = y0, times = times, func = ode_model, parms = c(Q = Q, PS = PS, fu = fu),
                      method = "lsoda")
  as.data.frame(out)
}

# ---- data ----------------------------------------------------------------------
clean_conc <- function(x) {
  if (is.character(x)) {
    x <- ifelse(grepl("<\\s*0?\\.?\\d+", x), 0, x)       # "< LLOQ" -> 0
    x <- suppressWarnings(as.numeric(x))
  }
  x
}

read_and_prepare_sheet <- function(file, sheet_name) {
  dat <- read_xlsx(file, sheet = sheet_name)
  dat$conc <- clean_conc(dat$conc)
  dat <- dat[!is.na(dat$time_h) & !is.na(dat$conc), ]
  dat <- dat[order(dat$time_h), ]
  meta <- data.frame(Group = unique(na.omit(as.character(dat$pdms_status))),
                     compound = unique(na.omit(as.character(dat$compound))),
                     phase = unique(na.omit(as.character(dat$phase))),
                     conc_uM = unique(na.omit(as.numeric(dat$conc_uM))))
  stopifnot(nrow(meta) == 1)
  dat_mean <- aggregate(conc ~ time_h, data = dat, FUN = mean)
  dat_mean <- dat_mean[order(dat_mean$time_h), ]
  uniq_t <- sort(unique(dat_mean$time_h))
  if (meta$phase %in% c("Dose1", "Dose2")) {
    # start the phase at the candidate time point with the highest mean concentration
    ok <- DROP_CANDIDATES[DROP_CANDIDATES + 1 <= length(uniq_t)]
    if (length(ok) == 0) stop("Too few timepoints for drop_candidates")
    cand <- do.call(rbind, lapply(ok, function(n_drop) {
      t_start <- uniq_t[n_drop + 1]
      data.frame(n_drop = n_drop, t_start = t_start, c0 = dat_mean$conc[match(t_start, dat_mean$time_h)])
    }))
    cand <- cand[order(-cand$c0, cand$n_drop), ]
    n_drop <- cand$n_drop[1]
    t_shift <- cand$t_start[1]
  } else if (meta$phase == "Washout") {
    if (DROP_WASHOUT_N + 1 > length(uniq_t)) stop("Too few washout timepoints")
    n_drop <- DROP_WASHOUT_N
    t_shift <- uniq_t[DROP_WASHOUT_N + 1]
  } else {
    stop("Unknown phase in sheet name/data")
  }
  dat2 <- dat[dat$time_h >= t_shift, ]
  dat_mean2 <- dat_mean[dat_mean$time_h >= t_shift, ]
  dat2$time_fit <- dat2$time_h - t_shift
  dat_mean2$time_fit <- dat_mean2$time_h - t_shift
  list(meta = meta, raw = dat2, mean = dat_mean2, c0 = dat_mean2$conc[dat_mean2$time_fit == 0][1],
       n_drop = n_drop, t_shift = t_shift, file = basename(file), sheet_name = sheet_name)
}

conditions <- list()
for (xlsx_file in XLSX_FILES) {
  sheets <- setdiff(excel_sheets(xlsx_file), c("README", "Readme", "readme"))
  sheets <- sheets[!(sheets %in% EXCLUDE_SHEETS)]
  for (sheet_name in sheets) {
    obj <- read_and_prepare_sheet(xlsx_file, sheet_name)
    if (!is.na(obj$meta$conc_uM) && obj$meta$conc_uM %in% SKIP_CONC_UM) next
    key <- paste(obj$meta$Group, obj$meta$compound, obj$meta$conc_uM, sep = "__")
    if (is.null(conditions[[key]])) conditions[[key]] <- list(meta = obj$meta, Dose1 = NULL, Dose2 = NULL, Washout = NULL)
    conditions[[key]][[obj$meta$phase]] <- obj
  }
}

# ---- fit -------------------------------------------------------------------------
fit_one_condition <- function(cond, fit_mode, fit_k3) {
  has_dose1 <- !is.null(cond$Dose1)
  has_dose2 <- !is.null(cond$Dose2)
  has_wash  <- !is.null(cond$Washout)
  if (!(has_dose1 || has_dose2)) stop("Condition has no Dose1 or Dose2 data")

  par_init <- PAR_INIT; lower <- PAR_LOWER; upper <- PAR_UPPER
  if (has_dose1 && fit_k3) {
    par_init <- c(par_init, k3_d1 = K3_INIT); lower <- c(lower, k3_d1 = K3_LOWER); upper <- c(upper, k3_d1 = K3_UPPER)
  }
  if (has_dose2 && fit_k3) {
    par_init <- c(par_init, k3_d2 = K3_INIT); lower <- c(lower, k3_d2 = K3_LOWER); upper <- c(upper, k3_d2 = K3_UPPER)
  }
  get_k3 <- function(par, phase) {
    if (!fit_k3) return(0)
    as.numeric(par[if (phase == "Dose1") "k3_d1" else "k3_d2"])
  }

  simulate_all_phases <- function(par) {
    fu <- as.numeric(par["fu"]); PS <- as.numeric(par["PS"])
    sims <- list()
    for (ph in c("Dose1", "Dose2")) {
      if (is.null(cond[[ph]])) next
      c0 <- cond[[ph]]$c0
      c3_0 <- get_k3(par, ph) * c0
      sims[[ph]] <- list(sim = simulate_phase(cond[[ph]]$mean$time_fit, c(c1 = c0, c2 = c0, c3 = c3_0), fu, PS),
                         c3_init = c3_0)
    }
    if (has_wash) {
      prev <- if (!is.null(sims$Dose2)) sims$Dose2$sim else sims$Dose1$sim
      c3_0w <- prev[nrow(prev), "c3"]
      c12_0w <- if (WASHOUT_C12_INIT == "zero") 0 else cond$Washout$c0
      sims$Washout <- list(sim = simulate_phase(cond$Washout$mean$time_fit, c(c1 = c12_0w, c2 = c12_0w, c3 = c3_0w), fu, PS),
                           c3_init = c3_0w)
    }
    sims
  }

  resfun <- function(par) {
    sims <- simulate_all_phases(par)
    res <- c()
    if (has_dose1) res <- c(res, cond$Dose1$mean$conc - sims$Dose1$sim$c1)
    if (has_dose2) res <- c(res, cond$Dose2$mean$conc - sims$Dose2$sim$c1)
    if (fit_mode == "joint" && has_wash) res <- c(res, cond$Washout$mean$conc - sims$Washout$sim$c1)
    res
  }

  fit <- minpack.lm::nls.lm(par = par_init, fn = resfun, lower = lower, upper = upper,
                            control = minpack.lm::nls.lm.control(maxiter = MAXITER))
  cf <- coef(fit)
  sims <- simulate_all_phases(cf)

  pred <- list(); obs_mean <- list(); obs_raw <- list()
  for (ph in c("Dose1", "Dose2", "Washout")) {
    if (is.null(cond[[ph]]) || is.null(sims[[ph]])) next
    used <- !(ph == "Washout" && fit_mode == "validate_washout")
    pred[[ph]] <- data.frame(sims[[ph]]$sim, phase = ph, time_fit = sims[[ph]]$sim$time, included_in_fit = used)
    obs_mean[[ph]] <- data.frame(cond[[ph]]$mean, phase = ph, included_in_fit = used)
    obs_raw[[ph]] <- data.frame(time_fit = cond[[ph]]$raw$time_fit, conc = cond[[ph]]$raw$conc,
                                phase = ph, included_in_fit = used)
  }
  pred <- bind_rows(pred); obs_mean <- bind_rows(obs_mean); obs_raw <- bind_rows(obs_raw)

  pvo <- left_join(obs_mean, pred[, c("phase", "time_fit", "c1")], by = c("phase", "time_fit"))
  rss_fit <- sum((pvo$conc[pvo$included_in_fit] - pvo$c1[pvo$included_in_fit])^2, na.rm = TRUE)
  wash_rmse <- if (has_wash) {
    w <- pvo[pvo$phase == "Washout", ]
    sqrt(mean((w$conc - w$c1)^2, na.rm = TRUE))
  } else NA_real_

  row <- data.frame(Group = cond$meta$Group, compound = cond$meta$compound, conc_uM = cond$meta$conc_uM,
                    fu = as.numeric(cf["fu"]), PS = as.numeric(cf["PS"]),
                    k3_d1 = if ("k3_d1" %in% names(cf)) as.numeric(cf["k3_d1"]) else 0,
                    k3_d2 = if ("k3_d2" %in% names(cf)) as.numeric(cf["k3_d2"]) else 0,
                    c3_init_d1 = if (!is.null(sims$Dose1)) sims$Dose1$c3_init else NA_real_,
                    c3_init_d2 = if (!is.null(sims$Dose2)) sims$Dose2$c3_init else NA_real_,
                    c3_init_washout = if (!is.null(sims$Washout)) sims$Washout$c3_init else NA_real_,
                    rss_fit = rss_fit, washout_rmse = wash_rmse, fit_mode = fit_mode,
                    has_Dose1 = has_dose1, has_Dose2 = has_dose2, has_Washout = has_wash,
                    t_shift_Dose1 = if (has_dose1) cond$Dose1$t_shift else NA_real_,
                    t_shift_Dose2 = if (has_dose2) cond$Dose2$t_shift else NA_real_,
                    t_shift_Washout = if (has_wash) cond$Washout$t_shift else NA_real_,
                    n_drop_Dose1 = if (has_dose1) cond$Dose1$n_drop else NA_real_,
                    n_drop_Dose2 = if (has_dose2) cond$Dose2$n_drop else NA_real_,
                    n_drop_Washout = if (has_wash) cond$Washout$n_drop else NA_real_)
  list(row = row, pred = pred, obs_mean = obs_mean, obs_raw = obs_raw)
}

panel <- NULL
for (fit_name in names(FITS)) {
  cfg <- FITS[[fit_name]]
  pc_log("empty-chip NSB fits: ", fit_name)
  rows <- list(); preds <- list(); obs_means <- list(); obs_raws <- list()
  failures <- data.frame(condition_key = character(), message = character())
  for (key in names(conditions)) {
    out <- tryCatch(fit_one_condition(conditions[[key]], cfg$mode, cfg$fit_k3), error = function(e) e)
    if (inherits(out, "error")) {
      failures <- rbind(failures, data.frame(condition_key = key, message = conditionMessage(out)))
      next
    }
    rows[[key]] <- out$row
    preds[[key]] <- cbind(out$pred, condition_key = key)
    obs_means[[key]] <- cbind(out$obs_mean, condition_key = key)
    obs_raws[[key]] <- cbind(out$obs_raw, condition_key = key)
    if (fit_name == PANEL_FIT && key == PANEL_CONDITION) panel <- out
  }
  write.csv(bind_rows(rows), pc_out(paste0("fig6d_fit_results_", fit_name, ".csv")), row.names = FALSE)
  write.csv(bind_rows(preds), pc_out(paste0("fig6d_predictions_", fit_name, ".csv")), row.names = FALSE)
  write.csv(bind_rows(obs_means), pc_out(paste0("fig6d_observed_mean_", fit_name, ".csv")), row.names = FALSE)
  write.csv(bind_rows(obs_raws), pc_out(paste0("fig6d_observed_raw_", fit_name, ".csv")), row.names = FALSE)
  write.csv(failures, pc_out(paste0("fig6d_fit_failures_", fit_name, ".csv")), row.names = FALSE)
  r <- bind_rows(rows)
  r <- r[r$Group == "PDMS" & r$compound == "DCF" & r$conc_uM == 100, ]
  pc_log(sprintf("  PDMS DCF 100 µM: K_media/plastic = %.4f, PS_NSB = %.5f mL/h", r$fu, r$PS))
}

# data of the Fig. 6D panel (PANEL_CONDITION, PANEL_FIT)
write.csv(panel$pred, pc_out("fig6d_panel_predictions.csv"), row.names = FALSE)
write.csv(panel$obs_mean, pc_out("fig6d_panel_observed_mean.csv"), row.names = FALSE)
write.csv(panel$obs_raw, pc_out("fig6d_panel_observed_raw.csv"), row.names = FALSE)
save_panel(plot_fig6d(panel$pred, panel$obs_mean, panel$obs_raw, base_size = 12), "fig6d_panel", width = 8, height = 3)
pc_log("done: ", pc_run_dir())
