# =============================================================================
# R/fiau_tktd.R - functions of the SI S7 analysis (FIAU toxicokinetics and
# toxicodynamics, CN Bio PhysioMimix LC12 liver plates)
#
# Used by scripts/si_s7_fiau_tktd.R.
#   data      fiau_read_lcms(), fiau_lcms_source(), fiau_read_albumin(),
#             fiau_vehicle_table(), fiau_dose_means()
#   exposure  FIAU_EVENTS, fiau_segments(), fiau_pc(), fiau_conc(), fiau_auc(),
#             fiau_cmax(), fiau_fit_pk()
#   TKTD      fiau_damage_at(), fiau_predictor(), fiau_fit_tktd(),
#             fiau_profile_ci(), fiau_rse(), fiau_fit_stats()
#   static    fiau_hill_fit(), fiau_hill_ci()
#
# Time axis: t = 0 h is the start of the FIAU treatment, 27.03.2023 (protocol day
# d4, Excel date serial 45012); the last albumin sample is day 14, 240 h.
# "Plate" = "chip" of the protocol: one LC12 plate with 12 liver wells; two wells
# per dose. Plate numbers follow the protocol (Chip 1 = lot 0027-007, Chip 2 =
# 0036-007, Chip 3 = 0037-005).
# =============================================================================

suppressPackageStartupMessages({ library(dplyr); library(minpack.lm) })

FIAU_LLOQ_UM   <- 0.005     # LLOQ of the LC-MS assay (µM)
FIAU_T0_SERIAL <- 45012     # Excel date serial of 27.03.2023 = treatment start
FIAU_TEND_H    <- 240       # end of the experiment (06.04.2023, d14)

# Dosing and exchange schedule (study protocol).
# type "full": the medium of the well is exchanged completely (1.8 mL);
# type "half": 900 µl are taken and 990 µl of fresh medium are added (plate 3, d12-d13);
# dose = multiple of the nominal dose in the new medium (0 = FIAU-free medium).
#   plate 1: 48-h treatments from d6 (29.03.) and d11 (03.04.), FIAU-free exchange on
#            d8 (31.03.), last treatment from d13 (05.04.) until the end; no FIAU before d6
#   plate 2: 8-h pulses on d4, d6, d8, d11, d13 (each followed by a FIAU-free exchange)
#   plate 3: continuous from d4; re-dosed (full exchange) on d7 (72 h); FIAU-free from
#            d8 (96 h); daily dosing from d11 (full exchange on d11, half exchanges on d12, d13)
FIAU_EVENTS <- list(
  `1` = data.frame(t = c(48, 96, 168, 216), type = "full", dose = c(1, 0, 1, 1)),
  `2` = data.frame(t = c(0, 8, 48, 56, 96, 104, 168, 176, 216, 224), type = "full",
                   dose = c(1, 0, 1, 0, 1, 0, 1, 0, 1, 0)),
  `3` = data.frame(t = c(0, 72, 96, 168, 192, 216), type = c("full", "full", "full", "full", "half", "half"),
                   dose = c(1, 1, 0, 1, 1, 1))
)
# nominal dose sets per plate (µM); plate 3 also had a cell-free well dosed with 100 µM
FIAU_DOSES <- list(`1` = c(0.2, 1, 10, 50, 100), `2` = c(1, 5, 10, 50, 250), `3` = c(0.1, 1, 10, 50))

# ---- LC-MS data ---------------------------------------------------------------------
# Analytical blocks ("Platte1".."Platte4") of the LC-MS report -> plates: the blocks
# are identified by their dose sets (checked below) and sampling dates.
FIAU_BLOCK_TO_PLATE <- c(Platte1 = 3L, Platte2 = 3L, Platte3 = 1L, Platte4 = 2L)

# Sampling dates of the sample names -> hours after the start of treatment
FIAU_DATE_TO_H <- c("27.03.2023 (d4)" = 0, "28.03.2023 (d5)" = 24, "29.03.2023 (d6)" = 48,
                    "30.03.2023 (d7)" = 72, "31.03.2023 (d8)" = 96, "03.04.2023 (d11)" = 168,
                    "04.04.2023 (d12)" = 192, "05.04.2023 (d13)" = 216, "06.04.2023 (d14)" = 240)

# LC-MS concentrations (data/si_s7_fiau/fiau_lc12_lcms_concentrations.csv: one row per
# sample, in the order of the LC-MS report; see data/README.md) -> list(samples =
# culture-medium samples, dosing = dosing-solution QC samples of 27.03.2023).
# Concentrations in µM; "< 0.005 [LLOQ]" -> blq.
fiau_read_lcms <- function(csv, lloq = FIAU_LLOQ_UM) {
  d <- utils::read.csv(csv, colClasses = "character", na.strings = "", check.names = FALSE,
                       fileEncoding = "UTF-8")
  d$block <- d$analytical_block
  parts <- strsplit(d$sample_name, ",")
  fld <- function(i) vapply(parts, function(x) if (length(x) >= i) trimws(x[[i]]) else NA_character_, "")
  d$vial <- trimws(d$vial)
  d$dose_str <- fld(1); d$date <- fld(3); d$tlab <- fld(4)
  d$vehicle <- grepl("DMSO", d$dose_str, ignore.case = TRUE)
  d$nominal_uM <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", d$dose_str)))
  d$nominal_uM[d$vehicle] <- 0
  txt <- trimws(d$conc_report)
  d$blq <- grepl("<", txt)
  d$conc_uM <- suppressWarnings(as.numeric(txt))               # NA for "< LLOQ"
  stopifnot(!any(is.na(d$conc_uM) & !d$blq), all(d$block %in% names(FIAU_BLOCK_TO_PLATE)))
  d$plate <- unname(FIAU_BLOCK_TO_PLATE[d$block])
  d$source <- fiau_lcms_source(d)                             # calibration range that gave the value

  is_dosing <- d$date == "45012"                              # dosing-solution QC (t = 0)
  dosing <- d[is_dosing, ]
  d <- d[!is_dosing, ]
  tl <- suppressWarnings(as.numeric(d$tlab))                  # plate 3, first 72 h: hours in the name
  d$time_h <- ifelse(!is.na(tl), tl, unname(FIAU_DATE_TO_H[d$date]) + ifelse(grepl("\\+8", d$tlab), 8, 0))
  stopifnot(!any(is.na(d$time_h)))
  # the blocks must carry the dose sets of the protocol (plate 3 + the cell-free 100 µM well)
  for (pl in 1:3) {
    got <- sort(unique(d$nominal_uM[d$plate == pl & !d$vehicle]))
    want <- sort(c(FIAU_DOSES[[as.character(pl)]], if (pl == 3) 100))
    stopifnot(isTRUE(all.equal(got, want)))
  }
  d <- d[order(d$plate, d$nominal_uM, d$time_h), ]            # stable: replicates keep the sheet order
  d <- d %>% dplyr::group_by(plate, nominal_uM, time_h) %>% dplyr::mutate(replicate = dplyr::row_number()) %>% dplyr::ungroup()
  d$well_id <- paste0("P", d$plate, "_", d$nominal_uM, "uM_r", d$replicate)
  d$cell_free <- d$plate == 3 & d$nominal_uM == 100          # liver well screwed out (protocol p. 2)
  # quality flag fixed before any fit: treated sample above 3 x nominal or vehicle above LLOQ
  d$flag_contam <- (!d$vehicle & !is.na(d$conc_uM) & d$conc_uM > 3 * d$nominal_uM) | (d$vehicle & !d$blq)
  d$conc_lloq2_uM <- ifelse(d$blq, lloq / 2, d$conc_uM)
  list(samples = as.data.frame(d), dosing = as.data.frame(dosing))
}

# Which calibration range of the assay (high or low) gave each reported value? `s` holds
# the reported values (blq, conc_uM) and the results of both calibration ranges
# (conc_high_range, conc_low_range; text as in the LC-MS report).
fiau_lcms_source <- function(s) {
  hi <- suppressWarnings(as.numeric(s$conc_high_range))
  lo <- suppressWarnings(as.numeric(s$conc_low_range))
  eq <- function(a, b) !is.na(a) & !is.na(b) & abs(a - b) < 1e-9 * pmax(1, abs(a))
  src <- ifelse(s$blq, "BLQ", ifelse(eq(s$conc_uM, hi), "high",
                               ifelse(eq(s$conc_uM, lo), "low", "not reproduced")))
  src
}

# ---- albumin -------------------------------------------------------------------------
# Sheet "Albumin + ALT" of the albumin / ALT workbook -> one row per well and sampling time.
# Times come from the date serial (chip 3 was sampled on other days than chips 1 and 2);
# the label "+ 8h" marks the end-of-pulse samples of plate 2.
fiau_read_albumin <- function(xlsx, t0 = FIAU_T0_SERIAL) {
  a <- as.data.frame(readxl::read_excel(xlsx, sheet = "Albumin + ALT", col_names = FALSE, .name_repair = "minimal"))
  colnames(a) <- c("date", "label", "chip", "conc", "ALT", "Alb"); a <- a[-1, ]
  a$Alb <- suppressWarnings(as.numeric(as.character(a$Alb)))
  a$serial <- suppressWarnings(as.numeric(as.character(a$date)))
  a$label <- trimws(as.character(a$label))
  a$plus8 <- grepl("8h", a$label)
  a$time_h <- (a$serial - t0) * 24 + ifelse(a$plus8, 8, 0)
  a$plate <- as.integer(gsub("Chip ", "", trimws(a$chip)))
  a$vehicle <- grepl("DMSO|ctrl", a$conc, ignore.case = TRUE)
  a$cell_free <- grepl("empty", a$conc, ignore.case = TRUE)
  a$nominal_uM <- ifelse(a$vehicle, 0, suppressWarnings(as.numeric(gsub("[^0-9.]", "", gsub(",", ".", a$conc)))))
  a$nominal_uM[a$cell_free] <- 100
  stopifnot(!any(is.na(a$time_h)), !any(is.na(a$plate)), !any(is.na(a$nominal_uM)))
  a <- a %>% dplyr::group_by(plate, time_h, nominal_uM) %>% dplyr::mutate(replicate = dplyr::row_number()) %>% dplyr::ungroup()
  a$well_id <- paste0("P", a$plate, "_", a$nominal_uM, "uM_r", a$replicate)
  as.data.frame(a[, c("plate", "well_id", "replicate", "nominal_uM", "time_h", "label", "plus8", "vehicle",
                      "cell_free", "Alb", "ALT")])
}

# Vehicle (0.1% DMSO) denominator per plate and time with the pre-specified flags:
#   copied   : the vehicle mean is identical to that of another time of the same plate
#   collapse : the vehicle mean is below `frac` of the first vehicle mean of the same plate and
#              sample type (regular samples vs the "+ 8 h" samples of plate 2)
#   keep_robust = neither (secondary, robustness analysis)
fiau_vehicle_table <- function(alb, frac = 0.30) {
  v <- alb[alb$vehicle & !is.na(alb$Alb), ] %>%
    dplyr::group_by(plate, time_h, plus8) %>%
    dplyr::summarise(n_veh = dplyr::n(), alb_veh = mean(Alb), veh_values = paste(round(Alb, 1), collapse = "/"),
                     .groups = "drop") %>%
    dplyr::arrange(plate, plus8, time_h) %>% dplyr::group_by(plate) %>%
    dplyr::mutate(copied = duplicated(alb_veh) | duplicated(alb_veh, fromLast = TRUE)) %>%
    dplyr::group_by(plate, plus8) %>%
    dplyr::mutate(first_veh = dplyr::first(alb_veh), veh_fraction = alb_veh / first_veh,
                  collapse = veh_fraction < frac) %>% dplyr::ungroup() %>%
    dplyr::mutate(keep_robust = !copied & !collapse) %>% dplyr::arrange(plate, time_h)
  as.data.frame(v)
}

# Normalised albumin: every treated well divided by the time-matched vehicle mean of the
# same plate (cell-free wells excluded). Returns the well-level table.
fiau_normalise <- function(alb, veh) {
  w <- alb[!alb$vehicle & !alb$cell_free & !is.na(alb$Alb) & alb$nominal_uM > 0, ]
  w <- merge(w, veh[, c("plate", "time_h", "n_veh", "alb_veh", "copied", "collapse", "keep_robust")],
             by = c("plate", "time_h"), all.x = FALSE)
  w$alb_ratio <- w$Alb / w$alb_veh
  w[order(w$plate, w$nominal_uM, w$time_h, w$replicate), ]
}

fiau_dose_means <- function(wells) {
  m <- wells %>% dplyr::group_by(plate, nominal_uM, time_h) %>%
    dplyr::summarise(ratio_mean = mean(alb_ratio), ratio_sd = if (dplyr::n() > 1) sd(alb_ratio) else NA_real_,
                     n_wells = dplyr::n(), n_veh = dplyr::first(n_veh), vehicle_copied = dplyr::first(copied),
                     vehicle_collapsed = dplyr::first(collapse), keep_robust = dplyr::first(keep_robust),
                     .groups = "drop")
  m <- as.data.frame(m); names(m)[names(m) == "ratio_mean"] <- "alb_ratio"; m
}

# ---- exposure: one-compartment medium model ---------------------------------------------
# Between events dC/dt = -k C (k = 0: FIAU stable; k < 0: net concentration increase, e.g.
# evaporation). At a full exchange C+ = f_res C- + (1 - f_res) F D; at a half exchange
# C+ = (0.9 C- + 0.99 F D)/1.89 (900 µl removed from 1.8 mL, 990 µl added). F = measured /
# nominal of the dosing medium, f_res = fraction of the old medium left behind by a full
# exchange. Segment table: C(t) = C0 exp(-k (t - t0)) for t0 < t <= t1. A sample taken at
# an event time is taken BEFORE the event (value of the segment that ends there).
fiau_segments <- function(plate, dose, k = 0, F = 1, f_res = 0, tend = FIAU_TEND_H) {
  ev <- FIAU_EVENTS[[as.character(plate)]]
  segs <- list(); cur_t <- 0; cur_C <- 0
  for (i in seq_len(nrow(ev))) {
    te <- ev$t[i]
    if (te > cur_t) {
      segs[[length(segs) + 1]] <- c(cur_t, te, cur_C)
      cur_C <- cur_C * exp(-k * (te - cur_t)); cur_t <- te
    }
    dn <- ev$dose[i] * dose * F
    cur_C <- if (ev$type[i] == "half") (0.9 * cur_C + 0.99 * dn) / 1.89 else f_res * cur_C + (1 - f_res) * dn
  }
  if (tend > cur_t) segs[[length(segs) + 1]] <- c(cur_t, tend, cur_C)
  s <- as.data.frame(do.call(rbind, segs)); names(s) <- c("t0", "t1", "C0"); s
}

# Concentration at the times `tout` (pre-event convention)
fiau_conc <- function(seg, tout, k = 0) {
  i <- findInterval(tout, seg$t0, left.open = TRUE)
  i[i == 0] <- 1
  seg$C0[i] * exp(-k * (tout - seg$t0[i]))
}

# Piecewise-constant version of a segment table (for the exact damage integration); for
# k != 0 every segment is cut into pieces of at most `step` hours (midpoint value)
fiau_pc <- function(seg, k = 0, step = 1) {
  if (k == 0) return(data.frame(t0 = seg$t0, t1 = seg$t1, C = seg$C0))
  do.call(rbind, lapply(seq_len(nrow(seg)), function(i) {
    n <- max(1, ceiling((seg$t1[i] - seg$t0[i]) / step))
    e <- seq(seg$t0[i], seg$t1[i], length.out = n + 1)
    mid <- (e[-1] + e[-length(e)]) / 2
    data.frame(t0 = e[-length(e)], t1 = e[-1], C = seg$C0[i] * exp(-k * (mid - seg$t0[i])))
  }))
}

# Cumulative AUC (µM h) over [0, t] and running maximum of C over times s < t, for a
# piecewise-constant table (k = 0)
fiau_auc <- function(pc, t) vapply(t, function(x) sum(pc$C * pmax(0, pmin(pc$t1, x) - pmin(pc$t0, x))), 0)
fiau_cmax <- function(pc, t) vapply(t, function(x) { s <- pc$t0 < x; if (any(s)) max(pc$C[s]) else 0 }, 0)

# Fit of the exposure model to the measured medium concentrations (log scale). `dat`:
# plate, nominal_uM, time_h, conc_uM. Parameters: F, f_res per plate; k fixed (kfix) or free
# (signed, "k_net"). Returns the fit and the predictions.
fiau_pk_predict <- function(dat, k, F, f_res) {
  out <- numeric(nrow(dat)); key <- paste(dat$plate, dat$nominal_uM)
  for (u in unique(key)) {
    s <- which(key == u); pl <- dat$plate[s[1]]
    seg <- fiau_segments(pl, dat$nominal_uM[s[1]], k = k, F = F, f_res = f_res[[as.character(pl)]])
    out[s] <- fiau_conc(seg, dat$time_h[s], k = k)
  }
  out
}

fiau_fit_pk <- function(dat, kfix = 0, shared_fres = FALSE) {
  pl <- sort(unique(dat$plate))
  nm <- c(if (is.null(kfix)) "k", "lF", if (shared_fres) "q" else paste0("q", pl))
  unpack <- function(par) {
    k <- if (is.null(kfix)) par[["k"]] else kfix
    q <- if (shared_fres) rep(par[["q"]], length(pl)) else par[paste0("q", pl)]
    list(k = k, F = exp(par[["lF"]]), f_res = setNames(plogis(unname(q)), as.character(pl)))
  }
  rf <- function(par) { u <- unpack(par)
    log(pmax(fiau_pk_predict(dat, u$k, u$F, u$f_res), 1e-10)) - log(dat$conc_uM) }
  start <- setNames(c(if (is.null(kfix)) 0, log(1.05), rep(qlogis(0.08), length(nm) - 1 - is.null(kfix))), nm)
  lower <- setNames(c(if (is.null(kfix)) -0.02, -2, rep(-8, length(nm) - 1 - is.null(kfix))), nm)
  upper <- setNames(c(if (is.null(kfix)) 0.02, 2, rep(4, length(nm) - 1 - is.null(kfix))), nm)
  f <- minpack.lm::nls.lm(par = start, fn = rf, lower = lower, upper = upper,
                          control = minpack.lm::nls.lm.control(maxiter = 400, ftol = 1e-12, ptol = 1e-12))
  n <- nrow(dat); p <- length(nm); rss <- sum(f$fvec^2)
  # standard errors from a numerical Jacobian (parameters on the log / logit scale)
  r0 <- rf(f$par)
  J <- sapply(seq_along(f$par), function(j) { pp <- f$par; h <- 1e-5 * max(1, abs(pp[j])); pp[j] <- pp[j] + h; (rf(pp) - r0) / h })
  V <- try(rss / (n - p) * solve(crossprod(J)), silent = TRUE)
  se <- if (inherits(V, "try-error")) rep(NA_real_, p) else sqrt(diag(V))
  names(se) <- nm
  list(par = f$par, se = setNames(se, nm), rss = rss, n = n, p = p, unpack = unpack, resfun = rf,
       aic = n * log(rss / n) + 2 * p, sigma = sqrt(rss / (n - p)), kfix = kfix, shared_fres = shared_fres)
}

# ---- toxicodynamics: damage model ---------------------------------------------------------
#   dD/dt = a(C) (1 - D) - krep D,  a(C) = kdam C/(EC50 + C),  albumin ratio = 1 - Emax D
# For a piecewise-constant exposure the solution is exact in every segment:
#   D(t) = Dinf + (D0 - Dinf) exp(-(a + krep)(t - t0)),  Dinf = a/(a + krep).
# sat = FALSE: a(C) = kdam C (no saturation; the EC50 -> infinity limit with kdam/EC50 fixed)
fiau_damage_at <- function(pc, tout, kdam, krep, EC50, sat = TRUE) {
  ord <- order(tout); ti <- tout[ord]; nt <- length(ti); res <- numeric(nt); j <- 1; D0 <- 0
  for (r in seq_len(nrow(pc))) {
    t0 <- pc$t0[r]; t1 <- pc$t1[r]; C <- pc$C[r]
    a <- if (sat) kdam * C / (EC50 + C) else kdam * C
    lam <- a + krep; Dinf <- if (lam > 0) a / lam else 0
    while (j <= nt && ti[j] <= t1) { res[j] <- Dinf + (D0 - Dinf) * exp(-lam * (ti[j] - t0)); j <- j + 1 }
    D0 <- Dinf + (D0 - Dinf) * exp(-lam * (t1 - t0))
  }
  while (j <= nt) { res[j] <- D0; j <- j + 1 }
  out <- numeric(length(tout)); out[ord] <- res; out
}

# Prediction function par -> albumin ratio for the rows of `dat`; `pcs` = named list
# ("plate dose") of piecewise-constant exposures. avg = TRUE: albumin accumulates in the
# medium since the last medium change, so the sample shows the mean of 1 - Emax D over
# that window (structural sensitivity analysis).
fiau_predictor <- function(dat, pcs, avg = FALSE) {
  idx <- split(seq_len(nrow(dat)), paste(dat$plate, dat$nominal_uM))
  win0 <- function(pl, t) { te <- FIAU_EVENTS[[as.character(pl)]]$t; te <- te[te < t]; if (length(te)) max(te) else 0 }
  function(par) {
    kd <- exp(par[1]); kr <- exp(par[2]); ec <- exp(par[3]); em <- plogis(par[4])
    out <- numeric(nrow(dat))
    for (u in names(idx)) {
      s <- idx[[u]]; pc <- pcs[[u]]
      if (!avg) out[s] <- 1 - em * fiau_damage_at(pc, dat$time_h[s], kd, kr, ec)
      else out[s] <- vapply(s, function(i) {
        w0 <- win0(dat$plate[i], dat$time_h[i])
        if (dat$time_h[i] <= w0) 1 - em * fiau_damage_at(pc, dat$time_h[i], kd, kr, ec)
        else mean(1 - em * fiau_damage_at(pc, seq(w0, dat$time_h[i], length.out = 41), kd, kr, ec))
      }, 0)
    }
    out
  }
}

# bounds on (log kdam, log krep, log EC50, logit Emax); krep and Emax typically end at a bound
FIAU_LOWER <- c(log(1e-4), log(1e-9), log(1e-3), qlogis(0.05))
FIAU_UPPER <- c(log(1), log(1), log(1e4), qlogis(1 - 1e-6))
FIAU_STARTS <- expand.grid(lkdam = log(c(0.002, 0.01, 0.05)), lkrep = log(c(1e-5, 1e-3, 0.01)),
                           lEC50 = log(c(0.05, 1, 10, 100)), qEmax = qlogis(c(0.7, 0.95)))

# Multistart least squares (Levenberg-Marquardt); `dat` needs plate, nominal_uM, time_h, alb_ratio
fiau_fit_tktd <- function(dat, pcs, starts = FIAU_STARTS, avg = FALSE) {
  pr <- fiau_predictor(dat, pcs, avg)
  rf <- function(par) pr(par) - dat$alb_ratio
  best <- NULL
  for (i in seq_len(nrow(starts))) {
    f <- try(minpack.lm::nls.lm(par = unlist(starts[i, ]), fn = rf, lower = FIAU_LOWER, upper = FIAU_UPPER,
                                control = minpack.lm::nls.lm.control(maxiter = 500, ftol = 1e-12, ptol = 1e-12)),
             silent = TRUE)
    if (!inherits(f, "try-error") && (is.null(best) || sum(f$fvec^2) < sum(best$fvec^2) - 1e-12)) best <- f
  }
  par <- setNames(best$par, c("lkdam", "lkrep", "lEC50", "qEmax"))
  list(par = par, value = c(kdam = exp(par[[1]]), krep = exp(par[[2]]), EC50 = exp(par[[3]]), Emax = plogis(par[[4]])),
       rss = sum(best$fvec^2), n = nrow(dat), p = 4, resfun = rf, pred = pr(best$par), info = best$info)
}

# 95% profile-likelihood interval of parameter `idx` (F threshold, RSS-based): grid over the
# transformed parameter, the other three re-optimised
fiau_profile_ci <- function(fit, idx, lo, hi, ngrid = 45) {
  p <- fit$par; n <- fit$n; np <- fit$p
  thr <- fit$rss * (1 + qf(0.95, 1, n - np) / (n - np))
  gr <- seq(lo, hi, length.out = ngrid); pv <- rep(NA_real_, ngrid)
  for (i in seq_along(gr)) {
    gi <- gr[i]
    ff <- try(minpack.lm::nls.lm(par = p[-idx], lower = FIAU_LOWER[-idx], upper = FIAU_UPPER[-idx],
                                 fn = function(q) { full <- numeric(4); full[idx] <- gi; full[-idx] <- q; fit$resfun(full) },
                                 control = minpack.lm::nls.lm.control(maxiter = 300)), silent = TRUE)
    if (!inherits(ff, "try-error")) pv[i] <- sum(ff$fvec^2)
  }
  ok <- which(!is.na(pv) & pv <= thr)
  if (!length(ok)) return(list(lo = NA, hi = NA, open_lo = NA, open_hi = NA, grid = gr, rss = pv, thr = thr))
  i1 <- min(ok); i2 <- max(ok)
  xlo <- if (i1 > 1 && is.finite(pv[i1 - 1])) gr[i1 - 1] + (thr - pv[i1 - 1]) * (gr[i1] - gr[i1 - 1]) / (pv[i1] - pv[i1 - 1]) else gr[i1]
  xhi <- if (i2 < ngrid && is.finite(pv[i2 + 1])) gr[i2] + (thr - pv[i2]) * (gr[i2 + 1] - gr[i2]) / (pv[i2 + 1] - pv[i2]) else gr[i2]
  list(lo = xlo, hi = xhi, open_lo = i1 == 1, open_hi = i2 == ngrid, grid = gr, rss = pv, thr = thr)
}

# Profile intervals of the four parameters (natural scale) as a data frame
FIAU_PROFILE_RANGE <- list(c(log(1e-4), log(1)), c(log(1e-9), log(1)), c(log(1e-3), log(1e4)), c(qlogis(0.3), qlogis(1 - 1e-4)))
# A limit that falls on the edge of the profile grid is open: reported as NA, except the two
# natural limits (krep >= 0 below, Emax <= 1 above).
fiau_profile_all <- function(fit, ngrid = 45) {
  inv <- list(exp, exp, exp, plogis)
  do.call(rbind, lapply(1:4, function(i) {
    r <- fiau_profile_ci(fit, i, FIAU_PROFILE_RANGE[[i]][1], FIAU_PROFILE_RANGE[[i]][2], ngrid)
    lo <- inv[[i]](r$lo); hi <- inv[[i]](r$hi)
    if (isTRUE(r$open_lo)) lo <- if (i == 2) 0 else NA_real_
    if (isTRUE(r$open_hi)) hi <- if (i == 4) 1 else NA_real_
    data.frame(parameter = c("kdam", "krep", "EC50_dam", "Emax")[i], ci_lo = lo, ci_hi = hi,
               open_lo = r$open_lo, open_hi = r$open_hi)
  }))
}

# Standard errors / RSE (%) from the Jacobian at the optimum. Parameters at a bound (krep ->
# 0, Emax -> 1) make the full covariance singular, so they are fixed at the bound and the
# remaining parameters are re-fitted (reduced model); their RSE is NA.
fiau_rse <- function(fit) {
  v <- fit$value; p <- fit$par
  at_b <- c(FALSE, v[["krep"]] < 1e-6, FALSE, v[["Emax"]] > 0.999)
  pfix <- p
  if (at_b[2]) pfix[2] <- FIAU_LOWER[2]
  if (at_b[4]) pfix[4] <- FIAU_UPPER[4]
  free <- which(!at_b); rse <- rep(NA_real_, 4)
  rf <- function(q) { full <- pfix; full[free] <- q; fit$resfun(full) }
  f2 <- minpack.lm::nls.lm(par = p[free], fn = rf, lower = FIAU_LOWER[free], upper = FIAU_UPPER[free],
                           control = minpack.lm::nls.lm.control(maxiter = 500, ftol = 1e-12))
  q <- f2$par; r2 <- rf(q); rss2 <- sum(r2^2)
  J <- sapply(seq_along(q), function(j) { qq <- q; h <- 1e-5 * max(1, abs(q[j])); qq[j] <- qq[j] + h; (rf(qq) - r2) / h })
  V <- try(rss2 / (fit$n - length(q)) * solve(crossprod(J)), silent = TRUE)
  if (!inherits(V, "try-error")) {
    se <- sqrt(diag(V))
    for (j in seq_along(free)) rse[free[j]] <- if (free[j] == 4) 100 * se[j] * (1 - plogis(q[j])) else 100 * se[j]
  }
  data.frame(parameter = c("kdam", "krep", "EC50_dam", "Emax"), rse_pct = rse, at_bound = at_b)
}

# Fit-level statistics: RSS, residual SD (n - p), R2 about the mean, correlation of the
# residuals with log10(dose) and time
fiau_fit_stats <- function(fit, dat) {
  res <- dat$alb_ratio - fit$pred
  ct <- function(x) { z <- suppressWarnings(cor.test(x, res)); c(r = unname(z$estimate), p = z$p.value) }
  cd <- ct(log10(dat$nominal_uM)); cT <- ct(dat$time_h)
  data.frame(n_points = fit$n, n_par = fit$p, rss = fit$rss, resid_sd = sqrt(fit$rss / (fit$n - fit$p)),
             r2 = 1 - fit$rss / sum((dat$alb_ratio - mean(dat$alb_ratio))^2),
             r_resid_logdose = cd[["r"]], p_resid_logdose = cd[["p"]], r_resid_time = cT[["r"]], p_resid_time = cT[["p"]])
}

# ---- static Emax (Hill, Emax fixed at 1) per plate and time ---------------------------------
# effect = 1 - albumin ratio (clipped to 0..1); weights w (wells per dose level)
fiau_hill_fit <- function(dose, ratio, w = rep(1, length(dose))) {
  eff <- pmax(0, pmin(1, 1 - ratio)); ld <- log(dose)
  rf <- function(p) sqrt(w) * (plogis(exp(p[2]) * (ld - p[1])) - eff)
  lo <- c(log(1e-3), log(0.1)); hi <- c(log(1e5), log(8)); best <- NULL
  for (a in log(median(dose)) + log(c(0.1, 1, 10))) for (b in log(c(0.5, 2))) {
    f <- try(minpack.lm::nls.lm(par = c(a, b), fn = rf, lower = lo, upper = hi,
                                control = minpack.lm::nls.lm.control(maxiter = 200)), silent = TRUE)
    if (!inherits(f, "try-error") && (is.null(best) || sum(f$fvec^2) < sum(best$fvec^2) - 1e-12)) best <- f
  }
  if (is.null(best)) return(list(EC50 = NA_real_, hill = NA_real_, rss = NA_real_, n = length(dose),
                                 max_effect = max(eff), eff = eff, ld = ld, w = w))
  list(EC50 = exp(best$par[1]), hill = exp(best$par[2]), rss = sum(best$fvec^2), n = length(dose),
       max_effect = max(eff), eff = eff, ld = ld, w = w)
}

# 95% profile-likelihood interval of the apparent EC50 (F threshold with n - 2 df; the Hill
# coefficient is re-optimised); NA where the profile does not cross the threshold
fiau_hill_ci <- function(fit) {
  n <- fit$n
  if (n < 4 || !is.finite(fit$rss)) return(c(lo = NA_real_, hi = NA_real_))
  grid <- seq(log(1e-3), log(1e5), length.out = 800)
  rss <- vapply(grid, function(th) optimize(function(lh) sum(fit$w * (plogis(exp(lh) * (fit$ld - th)) - fit$eff)^2),
                                            c(log(0.1), log(8)))$objective, 0)
  thr <- min(c(rss, fit$rss)) * (1 + qf(0.95, 1, n - 2) / (n - 2))
  ok <- which(rss <= thr); i1 <- min(ok); i2 <- max(ok)
  c(lo = if (i1 > 1) exp(grid[i1]) else NA_real_, hi = if (i2 < length(grid)) exp(grid[i2]) else NA_real_)
}
