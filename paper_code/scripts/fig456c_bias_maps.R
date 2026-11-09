# =============================================================================
# fig456c_bias_maps.R - Figs. 4C, 5C, 6C: bias of CL_int and P when
# evaporation (4C), sampling (5C) or NSB (6C) is present in the data but
# omitted from the analysis model
#
# For each compound on a 10 x 10 grid of CL_int (0.03-0.20 mL/h) and PS
# (0.003-0.10 mL/h), log-spaced: simulate the liver profile with the process,
# sample it at 0, 1, 2, 4, 8, 24, 48, 72 h with 0.1% proportional noise, refit
# CL_int and PS with evaporation, sampling and NSB all removed, and report
# est/true (cases "Evap", "Samp", "NSB").
#
# Usage  : Rscript scripts/fig456c_bias_maps.R            (all three cases)
#          Rscript scripts/fig456c_bias_maps.R Evap NSB   (selected cases)
# Inputs : data/models/NSBpara2DAnalysis.rds
# Outputs: fig456c_bias_<case>.csv (PS, CL, ratio_CL, ratio_PS, reported units),
#          fig456c_<case>_ratio_CL.png, fig456c_<case>_ratio_PS.png
#          (Evap -> Fig. 4C, Samp -> Fig. 5C, NSB -> Fig. 6C)
# Runtime: ~4.5 min per case (100 fits)
# =============================================================================

.f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
.root <- Sys.getenv("PAPER_CODE_ROOT", if (length(.f)) dirname(dirname(.f)) else ".")
source(file.path(.root, "R", "setup.R")); pc_setup(.root)
suppressPackageStartupMessages(library(ggplot2))

# ---- parameters ----------------------------------------------------------------
MODEL_FILE <- pc_data("models", "NSBpara2DAnalysis.rds")
CASES <- c("Evap", "Samp", "NSB")
FINAL_TIME_H <- 120                      # simulated horizon (observations 0-72 h)
GRID_N <- 10
CL_RANGE <- c(0.03, 0.20)                # mL/h, log-spaced
PS_RANGE <- c(0.003, 0.10)               # mL/h, log-spaced
PROP_NOISE_SD <- 0.001                   # proportional noise of the virtual data
SEED <- 1                                # set.seed() before each grid
# processes present in the virtual data (one per case)
EVAP_KEV_APICAL <- 0.0025                # mL/h (Fig. 4C)
EVAP_KEV_LIVER  <- 0.005                 # mL/h (Fig. 4C)
SAMP_VOLUME     <- 0.05                  # mL (Fig. 5C)
SAMP_TIMES_LIVER  <- c(0, 1, 2, 4, 8, 24, 48, 72)  # last time not applied; t = 0 ignored
SAMP_TIMES_APICAL <- c(8, 24, 72)                  # -> withdrawals at 8 and 24 h
NSB_PS          <- 0.05                  # mL/h, basolateral PS_NSB (Fig. 6C); apical = 1/5
NSB_K           <- 1                     # K_media/plastic (Fig. 6C)

args <- commandArgs(trailingOnly = TRUE)
if (length(args)) CASES <- intersect(CASES, args)

# ---- model ---------------------------------------------------------------------
# Base model of the misspecification analyses (build_misspec_app(): horizon
# 72 h, observations 0-72 h), then all processes off, horizon 120 h, then the
# process of the case switched on.
bias_app <- function(case) {
  app <- build_misspec_app(MODEL_FILE, kev_apical = 0.0026, kev_liver = 0.005, ps_nsb = 0.05,
                           k_media_plastic = 1, nsb_ratio = 5, samp_vol = 0.05)
  app$CompData[[1]]$Kev <- 0
  app$CompData[[3]]$Kev <- 0
  app$CompData[[5]]$fu <- 1
  app$CompData[[3]]$fu <- 1
  app$PermData[[2]]$PS <- 0
  app$PermData[[3]]$PS <- app$PermData[[2]]$PS / app$nsbRatio
  app$CompData[[1]]$ExMedia[, 2] <- 0
  app$CompData[[3]]$ExMedia[, 2] <- 0
  app$FinalTime <- FINAL_TIME_H
  app$CompData[[1]]$DoseSingData <- app$CompData[[1]]$Vol * 1 * app$MolarMass

  if (case == "Evap") {
    app$CompData[[1]]$Kev <- EVAP_KEV_APICAL
    app$CompData[[3]]$Kev <- EVAP_KEV_LIVER
  } else if (case == "Samp") {
    app$CompData[[3]]$ExMedia <- exchange_table(SAMP_TIMES_LIVER, SAMP_VOLUME)
    app$CompData[[1]]$ExMedia <- exchange_table(SAMP_TIMES_APICAL, SAMP_VOLUME)
  } else if (case == "NSB") {
    app$CompData[[5]]$fu <- NSB_K
    app$CompData[[6]]$fu <- NSB_K
    app$PermData[[2]]$PS <- NSB_PS
    app$PermData[[3]]$PS <- app$PermData[[2]]$PS / app$nsbRatio
  }
  app
}

run_case <- function(case) {
  app <- bias_app(case)
  PS_seq <- 10^seq(log10(PS_RANGE[1]), log10(PS_RANGE[2]), length.out = GRID_N)
  CL_seq <- 10^seq(log10(CL_RANGE[1]), log10(CL_RANGE[2]), length.out = GRID_N)
  grid <- expand.grid(PS = PS_seq, CL = CL_seq)
  set.seed(SEED)
  res <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
    app_pt <- copy_app(app)
    app_pt$PermData[[IDX_PS_GUT]]$PS <- grid$PS[i]
    app_pt$CompData[[IDX_LIVER]]$CL  <- grid$CL[i]
    app_truth <- simulate_truth(app_pt, prop_noise_sd = PROP_NOISE_SD)
    fit <- tryCatch(fit_neglecting_all(app_truth, grid$CL[i], grid$PS[i]),
                    error = function(e) c(est_CL = NA_real_, est_PS = NA_real_))
    if (i %% 25 == 0) pc_log("  ", case, " ", i, "/", nrow(grid))
    data.frame(PS = grid$PS[i], CL = grid$CL[i],
               ratio_CL = unname(fit["est_CL"]) / grid$CL[i],
               ratio_PS = unname(fit["est_PS"]) / grid$PS[i])
  }))
  res$CL_uL_min_1e6 <- cl_to_uLmin1e6(res$CL)
  res$P_1e6_cm_s    <- ps_to_P_1e6cms(res$PS)
  res
}

# ---- run -------------------------------------------------------------------------
for (case in CASES) {
  pc_log("bias map: ", case)
  res <- run_case(case)
  write.csv(res, pc_out(paste0("fig456c_bias_", case, ".csv")), row.names = FALSE)
  for (z in c("ratio_CL", "ratio_PS")) {
    save_panel(plot_bias_map(res, z, base_size = 14, bar_height_mm = 40),
               paste0("fig456c_", case, "_", z), width = 6.5, height = 4.8)
  }
  i <- which.min(abs(res$CL - CL_RANGE[1]) + abs(res$PS - PS_RANGE[1]))
  pc_log(sprintf("  %s: ratio_CL %.3f-%.3f, ratio_PS %.3f-%.3f; at CL %.2f / PS %.3f mL/h: ratio_CL = %.4f",
                 case, min(res$ratio_CL), max(res$ratio_CL), min(res$ratio_PS), max(res$ratio_PS),
                 res$CL[i], res$PS[i], res$ratio_CL[i]))
}
pc_log("done: ", pc_run_dir())
