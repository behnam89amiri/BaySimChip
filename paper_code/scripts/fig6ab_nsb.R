# =============================================================================
# fig6ab_nsb.R - Fig. 6A/B: effect of non-specific binding (NSB) to the device
#
#   A: gut apical and liver profiles of the generic compound for different
#      PS_NSB and K_media/plastic
#   B: liver exposure metrics over a 10 x 10 grid of PS_NSB and
#      K_media/plastic
#
# Model: NSBpara2DAnalysis.rds + apical solid compartment; one solid
# compartment on each side, common K_media/plastic, apical PS_NSB = 1/5 of the
# basolateral PS_NSB (contact area ratio).
#
# Inputs : data/models/NSBpara2DAnalysis.rds
# Outputs: fig6a_profiles.csv, fig6b_exposure_grid.csv,
#          fig6a_<compartment>.png, fig6b_<metric>.png
# Runtime: ~5 min
# =============================================================================

.f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
.root <- Sys.getenv("PAPER_CODE_ROOT", if (length(.f)) dirname(dirname(.f)) else ".")
source(file.path(.root, "R", "setup.R")); pc_setup(.root)
suppressPackageStartupMessages(library(ggplot2))

# ---- parameters ----------------------------------------------------------------
MODEL_FILE <- pc_data("models", "NSBpara2DAnalysis.rds")
NSB_RATIO <- 5                        # basolateral / apical PS_NSB (interfacial area ratio)
FINAL_TIME_A_H <- 100                 # Fig. 6A
COND_A <- data.frame(                 # Fig. 6A curves
  ps_nsb = c(0, 0.05, 0.1, 0.2, 0.05, 0.05),  # mL/h
  k      = c(1, 1, 1, 1, 0.1, 3)              # K_media/plastic
)
FINAL_TIME_B_H <- 300                 # Fig. 6B (longer horizon so that T1/2 is reached)
PS_NSB_GRID_B <- seq(0, 0.2, length.out = 10)  # mL/h
K_GRID_B      <- seq(0.2, 2, length.out = 10)

# ---- model -------------------------------------------------------------------------
nsb_app <- function(final_time) {
  app <- load_model(MODEL_FILE)
  app <- add_apical_nsb(app)
  app$FinalTime <- final_time
  app$CompData[[3]]$CL <- 0.1
  app$PermData[[1]]$PS <- 0.01
  for (i in 1:4) app$CompData[[i]]$Vol <- c(0.35, 0.8, 1.6, 0.2)[i]
  app$CompData[[1]]$DoseSingData <- app$CompData[[1]]$Vol * 1 * app$MolarMass
  app$CompData[[1]]$Kev <- 0
  app$CompData[[3]]$Kev <- 0
  ConstructDoseMediaEx(app)
}

# ---- Fig. 6A -------------------------------------------------------------------------
app <- nsb_app(FINAL_TIME_A_H)
prof <- do.call(rbind, lapply(seq_len(nrow(COND_A)), function(i) {
  app_i <- set_nsb(copy_app(app), COND_A$k[i], COND_A$ps_nsb[i], NSB_RATIO)
  cbind(run = i, ps_nsb = COND_A$ps_nsb[i], k = COND_A$k[i], simulate_profiles(app_i, 1:4))
}))
write.csv(prof, pc_out("fig6a_profiles.csv"), row.names = FALSE)

prof$cond <- factor(prof$ps_nsb, levels = c(0, 0.05, 0.1, 0.2))
prof$lty  <- factor(prof$k, levels = c(1, 0.1, 3))
ps_labels <- c(expression(PS[NSB]*" = 0 mL/h"), expression(PS[NSB]*" = 0.05 mL/h"),
               expression(PS[NSB]*" = 0.1 mL/h"), expression(PS[NSB]*" = 0.2 mL/h"))
k_labels <- c(expression(K[media/plastic]*" = 1"), expression(K[media/plastic]*" = 0.1"),
              expression(K[media/plastic]*" = 3"))
for (j in c(1, 3)) {
  p <- plot_profiles(prof, j, title = c("1" = "Gut apical", "3" = "Liver")[[as.character(j)]],
                     ylab = lab("conc_uM"), colours = setNames(HUE4, levels(prof$cond)),
                     colour_labels = ps_labels,
                     linetypes = c("1" = "solid", "0.1" = "dotted", "3" = "dotdash"),
                     linetype_labels = k_labels,
                     legend_inside = if (j == 1) c(0.98, 0.98) else NULL, base_size = 14)
  save_panel(p, paste0("fig6a_", c("1" = "gut_apical", "3" = "liver")[[as.character(j)]]), width = 7, height = 4)
}

# ---- Fig. 6B -------------------------------------------------------------------------
app <- nsb_app(FINAL_TIME_B_H)
grid <- expand.grid(K = K_GRID_B, PS_NSB = PS_NSB_GRID_B)[, c("PS_NSB", "K")]
res <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
  app_i <- set_nsb(copy_app(app), grid$K[i], grid$PS_NSB[i], NSB_RATIO)
  cbind(grid[i, ], liver_exposure(app_i))
}))
write.csv(res, pc_out("fig6b_exposure_grid.csv"), row.names = FALSE)

titles <- list(AUC = lab("AUC"), Cmax = lab("Cmax"), Tmax = lab("Tmax"), T12 = lab("T12"))
for (z in names(titles)) {
  p <- plot_heatmap(res, "PS_NSB", "K", z, titles[[z]], lab("ps_nsb"), lab("k_media"),
                    base_size = 14, bar_height_mm = 40)
  save_panel(p, paste0("fig6b_", z), width = 6, height = 4)
}
pc_log("done: ", pc_run_dir())
