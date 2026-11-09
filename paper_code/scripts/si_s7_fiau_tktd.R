# =============================================================================
# si_s7_fiau_tktd.R - SI Section S7 (Fig. S5): toxicokinetic-toxicodynamic model
# of FIAU hepatotoxicity in CN Bio PhysioMimix LC12 liver plates (3 plates
# ("chips") with different regimens, 12 wells each, n = 2 wells per dose,
# days 0-10 of treatment). Functions are in R/fiau_tktd.R.
#
# 1. PK data: tidy table of the LC-MS concentrations (for every value the
#    calibration range of the assay is recorded); the cell-free well (plate 3, 100 µM) and
#    samples flagged as contaminated are not used for fitting.
# 2. Exposure: one-compartment medium model that follows the protocol dosing and
#    exchange schedule exactly (R/fiau_tktd.R, FIAU_EVENTS); fitted to the measured
#    medium values (log scale): first-order loss k (net, signed), measured/nominal
#    factor F and the carry-over fraction f_res of a full exchange (per plate).
# 3. Albumin: every treated well divided by the time-matched mean of the 0.1% DMSO
#    vehicle wells of the same plate; dose-level means.
# 4. TKTD: dD/dt = kdam C/(EC50_dam + C) (1 - D) - krep D, albumin ratio = 1 - Emax D,
#    solved exactly between exposure changes; joint Levenberg-Marquardt fit (multi-
#    start; log / logit parameters) to the dose-level means of all plates, with
#    Jacobian RSE and 95% profile-likelihood intervals.
#    primary    : all dose-level means
#    robust     : drops ratios whose vehicle denominator is a copy or has fallen below
#                 30% of the first vehicle value of the plate (same sample type)
#    well_level : all individual wells (sensitivity)
#    plus structural sensitivity fits and fits of the plates separately.
# 5. Static Emax (model-free): Hill fit (Emax = 1) of the albumin effect vs nominal
#    dose per plate and time, apparent EC50 with 95% profile interval, drift factors.
#
# Inputs : data/si_s7_fiau/fiau_lc12_lcms_concentrations.csv (LC-MS concentrations)
#          data/si_s7_fiau/fiau_lc12_albumin_alt.xlsx (sheet "Albumin + ALT")
# Outputs: si_s7_pk_tidy.csv, si_s7_pk_parameters.csv, si_s7_albumin_dosemeans.csv,
#          si_s7_tktd_parameters.csv, si_s7_tktd_fit_summary.csv, si_s7_static_ec50.csv,
#          si_s7_ec50_drift.csv, si_s7_tktd.png/.pdf (A: albumin ratio and model per
#          plate, B: apparent EC50 vs time, C: damage state of plate 1)
# Runtime: ~3 min
# =============================================================================

.f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
.root <- Sys.getenv("PAPER_CODE_ROOT", if (length(.f)) dirname(dirname(.f)) else ".")
source(file.path(.root, "R", "setup.R")); pc_setup(.root)
source(pc_path("R", "fiau_tktd.R"))
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(minpack.lm)
})
t_script <- Sys.time()

# ---- parameters ----------------------------------------------------------------
LCMS_FILE <- pc_data("si_s7_fiau", "fiau_lc12_lcms_concentrations.csv")
ALB_FILE  <- pc_data("si_s7_fiau", "fiau_lc12_albumin_alt.xlsx")
ROBUST_VEHICLE_FRACTION <- 0.30        # robust fit: vehicle denominator < 30% of the first one is dropped
PROFILE_POINTS <- 45                   # grid points per profile-likelihood interval of the TKTD fit
K_RANGE <- c(-0.006, 0.006)            # net first-order loss rate k (1/h) of the exposure model: profile range
K_POINTS <- 49
MIN_EFFECT_INTERPOLATED <- 0.5         # apparent EC50 counts as interpolated if the top dose gives >= 50% effect,
MIN_DOSES_INTERPOLATED <- 4            #   at least 4 dose levels constrain the 2 Hill parameters, and it lies in the tested range
MAX_EC50_FACTOR <- 10                  # EC50 > 10 x the highest dose: not estimable (no effect at any dose)
GRID_H <- seq(0, FIAU_TEND_H, by = 1)  # time grid of the model curves in the figure

# ---- 1. PK data (LC-MS report) --------------------------------------------------------
lc <- fiau_read_lcms(LCMS_FILE)
pk <- lc$samples
pk$used_in_fit <- !pk$vehicle & !pk$flag_contam & !pk$blq & !pk$cell_free
fitdat <- pk[pk$used_in_fit, ]
src <- c(lc$samples$source, lc$dosing$source)
pc_log(sprintf("LC-MS data: %d sample rows (%d culture medium, %d dosing solutions); calibration range: %s",
               length(src), nrow(pk), nrow(lc$dosing), paste(names(table(src)), table(src), collapse = ", ")))
pc_log(sprintf("PK fit: %d of %d medium samples (excluded: %d vehicle, %d cell-free well, %d flagged, %d < LLOQ)",
               nrow(fitdat), nrow(pk), sum(pk$vehicle), sum(pk$cell_free), sum(pk$flag_contam & !pk$vehicle),
               sum(pk$blq & !pk$vehicle & !pk$flag_contam)))

# ---- 2. exposure model ----------------------------------------------------------------
pkf <- fiau_fit_pk(fitdat, kfix = 0)                          # final: k = 0, f_res per plate
u <- pkf$unpack(pkf$par)
F_HAT <- u$F; FRES <- u$f_res
ci_w <- function(est, se, inv) c(inv(est - 1.96 * se), inv(est + 1.96 * se))
ciF <- ci_w(pkf$par[["lF"]], pkf$se[["lF"]], exp)
ciq <- lapply(1:3, function(i) ci_w(pkf$par[[paste0("q", i)]], pkf$se[[paste0("q", i)]], plogis))
# alternative structures (AIC on the log-scale residuals)
pk_sh   <- fiau_fit_pk(fitdat, kfix = 0, shared_fres = TRUE)
pk_kfr  <- fiau_fit_pk(fitdat, kfix = NULL)
pk_ksh  <- fiau_fit_pk(fitdat, kfix = NULL, shared_fres = TRUE)
# profile of the net loss rate k in the final structure (F and f_res re-optimised)
kgrid <- seq(K_RANGE[1], K_RANGE[2], length.out = K_POINTS)
krss <- vapply(kgrid, function(k) fiau_fit_pk(fitdat, kfix = k)$rss, 0)
kthr <- min(krss) * (1 + qf(0.95, 1, pkf$n - pkf$p - 1) / (pkf$n - pkf$p - 1))
kin <- which(krss <= kthr)
cross <- function(i, j) kgrid[i] + (kthr - krss[i]) * (kgrid[j] - kgrid[i]) / (krss[j] - krss[i])
k_lo <- if (min(kin) > 1) cross(min(kin) - 1, min(kin)) else NA
k_hi <- if (max(kin) < K_POINTS) cross(max(kin) + 1, max(kin)) else NA
k_hat <- pk_kfr$par[["k"]]
# direct estimate from the continuous phase of plate 3 (no medium change from 0 to 72 h)
rich <- fitdat[fitdat$plate == 3 & fitdat$time_h <= 72, ]
lm3 <- lm(log(conc_uM) ~ time_h + factor(nominal_uM), data = rich)
slope3 <- coef(lm3)[["time_h"]]; slope3_ci <- confint(lm3)["time_h", ]
# diagnostics of the final model
fitdat$pred_uM <- fiau_pk_predict(fitdat, 0, F_HAT, FRES)
lr <- log(fitdat$pred_uM / fitdat$conc_uM)
trough <- fitdat$pred_uM < 0.25 * fitdat$nominal_uM
# cell-free well: not used for fitting; compared with the prediction of the cell-containing wells
cf <- pk[pk$cell_free & !pk$blq, ]
cf$pred_uM <- fiau_pk_predict(cf, 0, F_HAT, FRES)
cf_lr <- log(cf$conc_uM / cf$pred_uM)
cf_k <- nls.lm(par = c(k = 0), fn = function(p) log(fiau_pk_predict(cf, p[["k"]], F_HAT, FRES)) - log(cf$conc_uM),
                control = nls.lm.control(maxiter = 100))$par[["k"]]
# dosing solutions prepared on 27.03.2023 (LC-MS QC samples, t = 0)
ds <- lc$dosing[lc$dosing$nominal_uM > 0, ]
ds_ratio <- ds$conc_uM / ds$nominal_uM
pc_log(sprintf("exposure model: F = %.3f, f_res = %s, k_net = %.4f /h (95%% profile %.4f to %.4f)",
               F_HAT, paste(sprintf("%.3f", FRES), collapse = " / "), k_hat, k_lo, k_hi))

mfe <- function(x) exp(median(abs(x)))
pkp <- rbind(
  data.frame(section = "final model (k = 0)", item = "F, measured/nominal of the dosing medium", estimate = F_HAT,
             ci_lo = ciF[1], ci_hi = ciF[2], rse_pct = 100 * pkf$se[["lF"]], note = "Wald interval on the log scale"),
  do.call(rbind, lapply(1:3, function(i) data.frame(
    section = "final model (k = 0)", item = sprintf("f_res plate %d, fraction of the old medium left after a full exchange", i),
    estimate = FRES[[i]], ci_lo = ciq[[i]][1], ci_hi = ciq[[i]][2], rse_pct = 100 * pkf$se[[paste0("q", i)]],
    note = "Wald interval on the logit scale"))),
  data.frame(section = "loss rate", item = "k_net, net first-order loss (1/h; negative = concentration increase)",
             estimate = k_hat, ci_lo = k_lo, ci_hi = k_hi, rse_pct = NA,
             note = sprintf("k free in the final structure; 95%% profile interval; k = 0 is used (AIC %.1f vs %.1f); half-life if any loss >= %.0f h",
                            pkf$aic, pk_kfr$aic, log(2) / k_hi)),
  data.frame(section = "loss rate", item = "k from the continuous phase of plate 3 (1-72 h, no medium change; 1/h)",
             estimate = -slope3, ci_lo = -slope3_ci[[2]], ci_hi = -slope3_ci[[1]], rse_pct = NA,
             note = sprintf("concentration change over 72 h: %+.0f%% (evaporation); n = %d", 100 * (exp(slope3 * 72) - 1), nrow(rich))),
  data.frame(section = "structure comparison", item = c("k = 0, f_res shared", "k = 0, f_res per plate (used)",
                                                          "k free, f_res per plate", "k free, f_res shared"),
             estimate = c(pk_sh$aic, pkf$aic, pk_kfr$aic, pk_ksh$aic), ci_lo = NA, ci_hi = NA, rse_pct = NA,
             note = sprintf("AIC; RSS(log) %s", paste(sprintf("%.2f", c(pk_sh$rss, pkf$rss, pk_kfr$rss, pk_ksh$rss)), collapse = " / "))),
  data.frame(section = "diagnostics", item = c("observations used", "residual SD (log scale)", "median fold error",
                                                 "within 2-fold (%)", "within 1.5-fold (%)", "median fold error at trough samples (n)"),
             estimate = c(nrow(fitdat), pkf$sigma, mfe(lr), 100 * mean(abs(lr) < log(2)), 100 * mean(abs(lr) < log(1.5)), mfe(lr[trough])),
             ci_lo = NA, ci_hi = NA, rse_pct = NA,
             note = c("", "", "", "", "", sprintf("n = %d (predicted < 25%% of nominal)", sum(trough)))),
  data.frame(section = "dosing solutions", item = "geometric mean measured/nominal of the 8 dosing solutions (27.03.2023)",
             estimate = exp(mean(log(ds_ratio))), ci_lo = min(ds_ratio), ci_hi = max(ds_ratio), rse_pct = NA,
             note = "ci_lo / ci_hi = range; the medium after dosing shows F (above)"),
  data.frame(section = "cell-free well (not fitted)", item = c("geometric mean observed/predicted", "k fitted to the cell-free well alone (1/h)"),
             estimate = c(exp(mean(cf_lr)), cf_k), ci_lo = NA, ci_hi = NA, rse_pct = NA,
             note = c(sprintf("n = %d; plate-3 parameters, nominal 100 uM", nrow(cf)), "F and f_res fixed; not significant if near 0")))
write.csv(pkp, pc_out("si_s7_pk_parameters.csv"), row.names = FALSE)

# tidy PK table with the model prediction (non-vehicle samples)
pk$pred_uM <- NA_real_
nv <- !pk$vehicle
pk$pred_uM[nv] <- fiau_pk_predict(pk[nv, ], 0, F_HAT, FRES)
write.csv(pk %>% transmute(plate, well_id, replicate, nominal_uM, time_h, conc_uM, conc_lloq2_uM, blq,
                           vehicle, cell_free, flag_contam, used_in_fit, pred_uM, assay_sheet = source,
                           vial, sample_date = date, analytical_block = block),
          pc_out("si_s7_pk_tidy.csv"), row.names = FALSE)

# ---- 3. albumin -----------------------------------------------------------------------
alb <- fiau_read_albumin(ALB_FILE)
veh <- fiau_vehicle_table(alb, ROBUST_VEHICLE_FRACTION)
wells <- fiau_normalise(alb, veh)
dm <- fiau_dose_means(wells)
pc_log(sprintf("albumin: %d wells, %d dose-level means (%d kept by the robustness rule); vehicle wells per time: %s",
               nrow(wells), nrow(dm), sum(dm$keep_robust),
               paste(names(table(veh$n_veh)), table(veh$n_veh), sep = " x ", collapse = ", ")))
dm_k <- dm[dm$keep_robust, ]

# ---- 4. TKTD fits -----------------------------------------------------------------------
make_pcs <- function(F, fres, k = 0) {
  l <- list()
  for (pl in 1:3) for (dd in unique(dm$nominal_uM[dm$plate == pl]))
    l[[paste(pl, dd)]] <- fiau_pc(fiau_segments(pl, dd, k = k, F = F, f_res = fres[[as.character(pl)]]), k = k)
  l
}
pcs <- make_pcs(F_HAT, FRES)
pr_all <- fiau_predictor(dm, pcs)

fits <- list(); stats <- list(); rows <- list()
add_fit <- function(label, analysis, fit, dat, ci = NULL, note = "") {
  rse <- fiau_rse(fit); st <- fiau_fit_stats(fit, dat)
  v <- fit$value
  tab <- data.frame(parameter = c("kdam", "krep", "EC50_dam", "Emax"), estimate = unname(v[c("kdam", "krep", "EC50", "Emax")]),
                    unit = c("1/h", "1/h", "µM", "-"), rse_pct = rse$rse_pct, at_bound = rse$at_bound)
  tab$ci_lo <- NA_real_; tab$ci_hi <- NA_real_
  if (!is.null(ci)) { tab$ci_lo <- ci$ci_lo; tab$ci_hi <- ci$ci_hi }
  hl <- data.frame(parameter = c("damage_half_time", "repair_half_time"), estimate = log(2) / v[c("kdam", "krep")],
                   unit = "h", rse_pct = NA, at_bound = c(FALSE, rse$at_bound[2]), ci_lo = NA_real_, ci_hi = NA_real_)
  if (rse$at_bound[2]) hl$estimate[2] <- NA_real_              # krep at 0: no repair half-time, only a lower limit
  if (!is.null(ci)) {
    hl$ci_lo <- log(2) / c(ci$ci_hi[1], ci$ci_hi[2])
    hl$ci_hi <- log(2) / c(ci$ci_lo[1], ci$ci_lo[2])
    if (ci$ci_lo[2] <= 1e-8) hl$ci_hi[2] <- NA_real_           # open towards krep = 0
  }
  out <- cbind(fit = label, analysis = analysis, rbind(tab, hl), n_points = st$n_points, rss = st$rss,
               resid_sd = st$resid_sd, r2 = st$r2, note = note)
  rows[[length(rows) + 1]] <<- out
  stats[[length(stats) + 1]] <<- cbind(fit = label, analysis = analysis, st, F_vs_joint = NA_real_, p_vs_joint = NA_real_, note = note)
  fits[[label]] <<- fit
  invisible(fit)
}

t_fit <- Sys.time()
fit_p <- fiau_fit_tktd(dm, pcs)
ci_p <- fiau_profile_all(fit_p, PROFILE_POINTS)
add_fit("primary (dose-level means, all points)", "primary", fit_p, dm, ci_p)
pc_log(sprintf("primary fit: kdam %.4f /h, krep %.1e, EC50_dam %.2f µM, Emax %.3f; RSS %.3f (%.0f s)",
               fit_p$value[["kdam"]], fit_p$value[["krep"]], fit_p$value[["EC50"]], fit_p$value[["Emax"]], fit_p$rss,
               as.numeric(difftime(Sys.time(), t_fit, units = "secs"))))

pcs_k <- pcs[unique(paste(dm_k$plate, dm_k$nominal_uM))]
fit_r <- fiau_fit_tktd(dm_k, pcs)
ci_r <- fiau_profile_all(fit_r, PROFILE_POINTS)
add_fit("robust (dose-level means, vehicle rule)", "secondary", fit_r, dm_k, ci_r,
        sprintf("drops %d of %d dose means (vehicle denominator copied or < %.0f%% of the first)", nrow(dm) - nrow(dm_k), nrow(dm),
                100 * ROBUST_VEHICLE_FRACTION))

fit_w <- fiau_fit_tktd(wells, pcs)
ci_w2 <- fiau_profile_all(fit_w, PROFILE_POINTS)
add_fit("well level (all wells)", "sensitivity", fit_w, wells, ci_w2, "individual wells instead of dose-level means")

# one parameter set for the three regimens? plates fitted separately (primary data)
sep <- lapply(1:3, function(pl) { d <- dm[dm$plate == pl, ]; f <- fiau_fit_tktd(d, pcs); f })
for (pl in 1:3) add_fit(sprintf("plate %d alone (dose-level means)", pl), "sensitivity", sep[[pl]], dm[dm$plate == pl, ], NULL,
                        "plate fitted separately")
rss_sep <- sum(vapply(sep, function(f) f$rss, 0)); df_sep <- nrow(dm) - 12
Fst <- ((fit_p$rss - rss_sep) / 8) / (rss_sep / df_sep)
pF <- pf(Fst, 8, df_sep, lower.tail = FALSE)
stats[[length(stats) + 1]] <- data.frame(fit = "three plates separately (12 parameters) vs joint (4)", analysis = "sensitivity",
  n_points = nrow(dm), n_par = 12, rss = rss_sep, resid_sd = sqrt(rss_sep / df_sep), r2 = NA, r_resid_logdose = NA,
  p_resid_logdose = NA, r_resid_time = NA, p_resid_time = NA, F_vs_joint = Fst, p_vs_joint = pF,
  note = "F test, 8 and n - 12 degrees of freedom")
pc_log(sprintf("one parameter set vs plate-specific sets: F(8, %d) = %.2f, p = %.2f", df_sep, Fst, pF))

# structural sensitivity of the exposure model and of the albumin readout (point estimates)
# sensitivity fits start from the primary optimum plus every 9th point of the start grid
starts_sens <- rbind(as.data.frame(t(pmin(pmax(fit_p$par, FIAU_LOWER + 1e-6), FIAU_UPPER - 1e-6))),
                     FIAU_STARTS[seq(1, nrow(FIAU_STARTS), by = 9), ])
sens <- function(label, pcs_s, note, avg = FALSE) {
  f <- fiau_fit_tktd(dm, pcs_s, starts = starts_sens, avg = avg)
  add_fit(label, "sensitivity", f, dm, NULL, note)
}
sens("exposure without carry-over (f_res = 0)", make_pcs(F_HAT, c(`1` = 0, `2` = 0, `3` = 0)),
     "FIAU-free periods set to exactly 0 after each exchange")
fr_sh <- pk_sh$unpack(pk_sh$par)
sens("exposure with one carry-over fraction", make_pcs(fr_sh$F, fr_sh$f_res), sprintf("f_res = %.3f for all plates", fr_sh$f_res[[1]]))
kn <- pk_kfr$unpack(pk_kfr$par)
sens("exposure with net concentration drift k_net", make_pcs(kn$F, kn$f_res, kn$k),
     sprintf("k_net = %.4f /h (evaporation); not significant", kn$k))
sens("albumin averaged over the accumulation window", pcs, "albumin sample = mean of the model ratio since the last medium change", avg = TRUE)
pc_log(sprintf("TKTD fits done (%.0f s)", as.numeric(difftime(Sys.time(), t_fit, units = "secs"))))

tktd <- do.call(rbind, rows); rownames(tktd) <- NULL
write.csv(tktd, pc_out("si_s7_tktd_parameters.csv"), row.names = FALSE)
write.csv(do.call(rbind, stats), pc_out("si_s7_tktd_fit_summary.csv"), row.names = FALSE)

dm$pred_primary <- pr_all(fit_p$par)
dm$resid_primary <- dm$alb_ratio - dm$pred_primary
dm$pred_robust <- pr_all(fit_r$par)
write.csv(dm, pc_out("si_s7_albumin_dosemeans.csv"), row.names = FALSE)

# ---- 5. static Emax (Hill, Emax = 1) per plate and time ---------------------------------
first_dose <- vapply(1:3, function(pl) { ev <- FIAU_EVENTS[[as.character(pl)]]; min(ev$t[ev$dose > 0]) }, 0)
ec <- list()
for (pl in 1:3) for (tt in sort(unique(dm$time_h[dm$plate == pl]))) {
  m <- dm[dm$plate == pl & dm$time_h == tt, ]
  wl <- wells[wells$plate == pl & wells$time_h == tt, ]
  r <- data.frame(plate = pl, time_h = tt, n_doses = nrow(m), n_wells = nrow(wl), EC50_uM = NA_real_, EC50_lo = NA_real_,
                  EC50_hi = NA_real_, hill = NA_real_, max_effect = max(pmax(0, pmin(1, 1 - m$alb_ratio))),
                  exposed = tt > first_dose[pl], keep_robust = m$keep_robust[1], note = "")
  if (nrow(m) < 3) r$note <- sprintf("%d dose levels: no fit", nrow(m))
  if (nrow(m) >= 3) {
    h <- fiau_hill_fit(wl$nominal_uM, wl$alb_ratio)
    if (is.finite(h$EC50) && h$EC50 <= MAX_EC50_FACTOR * max(m$nominal_uM)) {
      ci <- fiau_hill_ci(h)
      r$EC50_uM <- h$EC50; r$hill <- h$hill; r$EC50_lo <- ci[["lo"]]; r$EC50_hi <- ci[["hi"]]
    } else r$note <- "EC50 not estimable (no effect at any dose)"
    if (nrow(m) < MIN_DOSES_INTERPOLATED && !nzchar(r$note)) r$note <- sprintf("%d dose levels: not constrained", nrow(m))
  }
  r$in_range <- is.finite(r$EC50_uM) && r$EC50_uM >= min(m$nominal_uM) && r$EC50_uM <= max(m$nominal_uM)
  r$interpolated <- r$exposed && r$in_range && r$max_effect >= MIN_EFFECT_INTERPOLATED && nrow(m) >= MIN_DOSES_INTERPOLATED
  ec[[length(ec) + 1]] <- r
}
ec <- do.call(rbind, ec)
write.csv(ec, pc_out("si_s7_static_ec50.csv"), row.names = FALSE)

drift_of <- function(d, label) {
  d <- d[is.finite(d$EC50_uM), ]
  do.call(rbind, lapply(sort(unique(d$plate)), function(pl) {
    s <- d[d$plate == pl, ]; s <- s[order(s$time_h), ]
    if (nrow(s) < 2) return(data.frame(variant = label, plate = pl, n_fits = nrow(s), t_first_h = NA, EC50_first = NA,
                                       t_last_h = NA, EC50_last = NA, drift_factor = NA))
    data.frame(variant = label, plate = pl, n_fits = nrow(s), t_first_h = s$time_h[1], EC50_first = s$EC50_uM[1],
               t_last_h = s$time_h[nrow(s)], EC50_last = s$EC50_uM[nrow(s)], drift_factor = s$EC50_uM[1] / s$EC50_uM[nrow(s)])
  }))
}
drift <- rbind(
  drift_of(ec, "all estimable fits (pre-dose and extrapolated included)"),
  drift_of(ec[which(ec$EC50_uM < 1000), ], "all fits with EC50 < 1000 uM"),
  drift_of(ec[ec$exposed, ], "fits after the first dose"),
  drift_of(ec[ec$interpolated, ], "interpolated fits only (>= 4 dose levels, top dose >= 50% effect)"),
  drift_of(ec[ec$interpolated & ec$keep_robust, ], "interpolated fits, robust points only"))
write.csv(drift, pc_out("si_s7_ec50_drift.csv"), row.names = FALSE)
print(drift, digits = 4, row.names = FALSE)

# ---- 6. figure (Fig. S5, print size) ------------------------------------------------------
BASE <- 8
th <- theme_bw(base_size = BASE) +
  theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey92"),
        legend.key.size = unit(0.30, "cm"), legend.margin = margin(0, 0, 0, 0),
        plot.tag = element_text(face = "bold", size = BASE + 1), axis.title = element_text(size = BASE),
        legend.title = element_text(size = BASE - 1), legend.text = element_text(size = BASE - 1),
        strip.text = element_text(size = BASE - 1))
dose_levels <- sort(unique(dm$nominal_uM))
cols <- setNames(colorRampPalette(c("#1A9850", "#FEE08B", "#D73027", "#7B3294"))(length(dose_levels)), as.character(dose_levels))
lab_plate <- c(`1` = "Plate 1: 48-h treatments", `2` = "Plate 2: 8-h pulses", `3` = "Plate 3: continuous, then daily")
v <- fit_p$value
traj <- do.call(rbind, lapply(names(pcs), function(key) {
  pl <- as.integer(strsplit(key, " ")[[1]][1]); dd <- as.numeric(strsplit(key, " ")[[1]][2])
  D <- fiau_damage_at(pcs[[key]], GRID_H, v[["kdam"]], v[["krep"]], v[["EC50"]])
  data.frame(plate = pl, nominal_uM = dd, time_h = GRID_H, D = D, alb_pred = 1 - v[["Emax"]] * D)
}))
dm$panel <- lab_plate[as.character(dm$plate)]; traj$panel <- lab_plate[as.character(traj$plate)]
pA <- ggplot() +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey60", linewidth = 0.3) +
  geom_line(data = traj, aes(time_h, alb_pred, colour = factor(nominal_uM)), linewidth = 0.5) +
  geom_point(data = dm, aes(time_h, alb_ratio, colour = factor(nominal_uM), shape = keep_robust), size = 1.1) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), guide = "none") +
  facet_wrap(~panel, nrow = 1) + scale_colour_manual(values = cols, name = "FIAU (µM)") +
  coord_cartesian(ylim = c(0, 1.6)) + labs(x = "Time (h)", y = "Albumin\n(fraction of vehicle)") + th
ec$kind <- ifelse(!ec$exposed, "before first dose", ifelse(ec$interpolated, "interpolated", "extrapolated or < 4 doses"))
ecp <- ec[is.finite(ec$EC50_uM), ]
pB <- ggplot(ecp, aes(time_h, EC50_uM, colour = factor(plate), shape = kind)) +
  geom_hline(yintercept = max(dm$nominal_uM), linetype = "dotted", colour = "grey40", linewidth = 0.3) +
  geom_linerange(data = ecp[is.finite(ecp$EC50_lo) & is.finite(ecp$EC50_hi) & ecp$kind == "interpolated", ],
                 aes(ymin = EC50_lo, ymax = EC50_hi), linewidth = 0.3, show.legend = FALSE) +
  geom_line(data = ecp[ecp$kind == "interpolated", ], aes(group = plate), linewidth = 0.5) +
  geom_point(size = 1.6) +
  scale_shape_manual(values = c(interpolated = 16, `extrapolated or < 4 doses` = 2, `before first dose` = 1), name = NULL) +
  scale_y_log10(labels = scales::label_number(big.mark = "", drop0trailing = TRUE)) +
  scale_colour_brewer(palette = "Set1", name = "Plate") +
  coord_cartesian(ylim = c(0.3, 1500)) + labs(x = "Time (h)", y = "Apparent EC50 (µM)") + th
pC <- ggplot(traj[traj$plate == 1, ], aes(time_h, D, colour = factor(nominal_uM))) +
  geom_line(linewidth = 0.5) + scale_colour_manual(values = cols, name = "FIAU (µM)") +
  coord_cartesian(ylim = c(0, 1)) + labs(x = "Time (h)", y = "Damage state D (plate 1)") + th
fig <- pA / (pB | pC) + plot_layout(heights = c(1, 0.95)) + plot_annotation(tag_levels = "A")
save_figure(fig, "si_s7_tktd", width = 6.7, height = 5.0)

pc_log(sprintf("done in %.1f min: %s", as.numeric(difftime(Sys.time(), t_script, units = "mins")), pc_run_dir()))
