# =============================================================================
# fig5ab_sampling.R - Fig. 5A/B: effect of sampling
#
#   A: gut apical and liver profiles of the generic compound for sample
#      volumes 0 / 0.05 / 0.10 / 0.15 mL
#   B: liver exposure metrics versus the number of apical (1-7) and liver
#      (1-32) samples of 0.05 mL, taken at equal intervals over 0-98 h
#
# Inputs : data/models/Kev2DAnalysis.rds
# Outputs: fig5a_profiles.csv, fig5b_exposure_grid.csv, fig5b_T12_change.csv,
#          fig5a_<compartment>.png, fig5b_<metric>.png
# Runtime: ~3 min
# =============================================================================

.f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
.root <- Sys.getenv("PAPER_CODE_ROOT", if (length(.f)) dirname(dirname(.f)) else ".")
source(file.path(.root, "R", "setup.R")); pc_setup(.root)
suppressPackageStartupMessages(library(ggplot2))

# ---- parameters ----------------------------------------------------------------
MODEL_FILE <- pc_data("models", "Kev2DAnalysis.rds")
FINAL_TIME_H <- 100                         # Figs. 5A and 5B
SAMPLE_VOLUMES_A <- c(0, 0.05, 0.10, 0.15)  # mL, Fig. 5A
# Fig. 5A schedule. As in the original analysis, the last listed time (72 h) is
# not applied, and the engine ignores the t = 0 withdrawal: effective
# withdrawals are liver 1, 2, 4, 8, 24, 48 h and apical 8, 24 h.
TIMES_LIVER_A  <- c(0, 1, 2, 4, 8, 24, 48, 72)
TIMES_APICAL_A <- c(8, 24, 72)
SAMPLE_VOLUME_B <- 0.05                     # mL, Fig. 5B
LAST_SAMPLE_FRACTION_B <- 0.98              # samples at round(seq(0, 0.98 x FinalTime, length.out = N), 2)

# ---- Fig. 5A -------------------------------------------------------------------
app <- setup_generic_compound(MODEL_FILE, final_time = FINAL_TIME_H)
app$CompData[[1]]$Kev <- 0
app$CompData[[3]]$Kev <- 0
prof <- do.call(rbind, lapply(SAMPLE_VOLUMES_A, function(v) {
  app_i <- copy_app(app)
  app_i$CompData[[3]]$ExMedia <- exchange_table(TIMES_LIVER_A, v)
  app_i$CompData[[1]]$ExMedia <- exchange_table(TIMES_APICAL_A, v)
  app_i <- ConstructDoseMediaEx(app_i)
  cbind(cond = if (v == 0) "0 mL" else sprintf("%.2f mL", v), simulate_profiles(app_i))
}))
write.csv(prof, pc_out("fig5a_profiles.csv"), row.names = FALSE)

prof$cond <- factor(prof$cond, levels = c("0 mL", "0.05 mL", "0.10 mL", "0.15 mL"))
cols <- setNames(HUE4, levels(prof$cond))
for (j in c(1, 3)) {
  p <- plot_profiles(prof, j, title = c("1" = "Gut apical", "3" = "Liver")[[as.character(j)]],
                     ylab = lab("conc_au"), colours = cols, colour_title = "Sample volume",
                     legend_inside = if (j == 1) c(0.98, 0.95) else NULL, base_size = 14)
  save_panel(p, paste0("fig5a_", c("1" = "gut_apical", "3" = "liver")[[as.character(j)]]), width = 6, height = 4)
}

# ---- Fig. 5B -------------------------------------------------------------------
app <- setup_generic_compound(MODEL_FILE, final_time = FINAL_TIME_H)
app$CompData[[1]]$Kev <- 0
app$CompData[[3]]$Kev <- 0
n_apical_grid <- seq(1, app$CompData[[1]]$Vol / SAMPLE_VOLUME_B, by = 1)   # 1..7 (volume limit)
n_liver_grid  <- seq(1, app$CompData[[3]]$Vol / SAMPLE_VOLUME_B, by = 1)   # 1..32
grid <- expand.grid(N_liver = n_liver_grid, N_apical = n_apical_grid)[, c("N_apical", "N_liver")]
res <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
  app_i <- copy_app(app)
  t_api <- round(seq(0, LAST_SAMPLE_FRACTION_B * app_i$FinalTime, length.out = grid$N_apical[i]), 2)[-1]
  t_liv <- round(seq(0, LAST_SAMPLE_FRACTION_B * app_i$FinalTime, length.out = grid$N_liver[i]), 2)[-1]
  app_i$CompData[[1]]$ExMedia <- cbind(t_api, rep(SAMPLE_VOLUME_B, length(t_api)), rep(0, length(t_api)))
  app_i$CompData[[3]]$ExMedia <- cbind(t_liv, rep(SAMPLE_VOLUME_B, length(t_liv)), rep(0, length(t_liv)))
  app_i <- ConstructDoseMediaEx(app_i)
  cbind(grid[i, ], liver_exposure(app_i))
}))
write.csv(res, pc_out("fig5b_exposure_grid.csv"), row.names = FALSE)

# T1/2 change of the densest schedule (7 apical x 32 liver samples) vs 1 x 1
# (main text: "reduced the profile-based T1/2 by approximately 30%")
t12_1 <- res$T12[res$N_apical == 1 & res$N_liver == 1]
t12_d <- res$T12[res$N_apical == max(n_apical_grid) & res$N_liver == max(n_liver_grid)]
t12 <- data.frame(T12_1x1_h = t12_1, T12_densest_h = t12_d, change_pct = 100 * (t12_d / t12_1 - 1))
write.csv(t12, pc_out("fig5b_T12_change.csv"), row.names = FALSE)
print(t12)

titles <- list(AUC = lab("AUC"), Cmax = lab("Cmax"), Tmax = lab("Tmax"), T12 = lab("T12"))
for (z in names(titles)) {
  p <- plot_heatmap(res, "N_apical", "N_liver", z, titles[[z]], lab("n_apical"), lab("n_liver"),
                    base_size = 14, bar_height_mm = 40)
  save_panel(p, paste0("fig5b_", z), width = 6, height = 4)
}
pc_log("done: ", pc_run_dir())
