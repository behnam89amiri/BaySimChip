# =============================================================================
# fig4ab_evaporation.R - Fig. 4A/B: effect of medium evaporation
#
#   A: gut apical and liver profiles of the generic compound with and without
#      evaporation
#   B: liver exposure metrics (AUC, Cmax, Tmax, T1/2) over a 5 x 5 grid of
#      apical and liver evaporation rates
#
# Inputs : data/models/Kev2DAnalysis.rds
# Outputs: fig4a_profiles.csv, fig4a_liver_cmax.csv, fig4b_exposure_grid.csv,
#          fig4a_<compartment>.png, fig4b_<metric>.png
# Runtime: ~1.5 min
# =============================================================================

.f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
.root <- Sys.getenv("PAPER_CODE_ROOT", if (length(.f)) dirname(dirname(.f)) else ".")
source(file.path(.root, "R", "setup.R")); pc_setup(.root)
suppressPackageStartupMessages(library(ggplot2))

# ---- parameters ----------------------------------------------------------------
MODEL_FILE <- pc_data("models", "Kev2DAnalysis.rds")
# generic compound (CL_int 0.1 mL/h, PS 0.01 mL/h, 1 µM apical) and volumes
# 0.35 / 0.8 / 1.6 / 0.2 mL: setup_generic_compound() in R/scenarios.R
FINAL_TIME_A_H  <- 110     # Fig. 4A horizon
KEV_APICAL_A    <- 0.003   # mL/h, Fig. 4A (upper limit: measured medium volume loss)
KEV_LIVER_A     <- 0.0055  # mL/h, Fig. 4A (upper limit: measured medium volume loss)
FINAL_TIME_B_H  <- 120     # Fig. 4B horizon
N_GRID_B        <- 5       # Kev grid 0 ... V/FINAL_TIME_B_H (volume left > 0), 5 x 5

# ---- Fig. 4A -------------------------------------------------------------------
app <- setup_generic_compound(MODEL_FILE, final_time = FINAL_TIME_A_H)
app <- ConstructDoseMediaEx(app)
conds <- data.frame(cond = c("w/ evaporation", "w/o evaporation"),
                    kev_apical = c(KEV_APICAL_A, 0),
                    kev_liver = c(KEV_LIVER_A, 0))
prof <- do.call(rbind, lapply(seq_len(nrow(conds)), function(i) {
  app_i <- copy_app(app)
  app_i$CompData[[1]]$Kev <- conds$kev_apical[i]
  app_i$CompData[[3]]$Kev <- conds$kev_liver[i]
  cbind(cond = conds$cond[i], simulate_profiles(app_i))
}))
write.csv(prof, pc_out("fig4a_profiles.csv"), row.names = FALSE)

liver <- prof[prof$compartment == 3, ]
cmax <- do.call(rbind, lapply(split(liver, liver$cond), function(d) {
  data.frame(cond = d$cond[1], Cmax = max(d$conc), Tmax_h = d$time_h[which.max(d$conc)])
}))
write.csv(cmax, pc_out("fig4a_liver_cmax.csv"), row.names = FALSE)
print(cmax, digits = 6)

prof$cond <- factor(prof$cond, levels = c("w/o evaporation", "w/ evaporation"))
cols <- setNames(HUE2, levels(prof$cond))
for (j in c(1, 3)) {
  p <- plot_profiles(prof, j, title = c("1" = "Gut apical", "3" = "Liver")[[as.character(j)]],
                     ylab = lab("conc_au"), colours = cols,
                     legend_inside = if (j == 1) c(0.98, 0.95) else NULL, base_size = 14)
  save_panel(p, paste0("fig4a_", c("1" = "gut_apical", "3" = "liver")[[as.character(j)]]), width = 6.7, height = 4)
}

# ---- Fig. 4B -------------------------------------------------------------------
app <- setup_generic_compound(MODEL_FILE, final_time = FINAL_TIME_B_H)
kev1_grid <- seq(0, app$CompData[[1]]$Vol / app$FinalTime, length.out = N_GRID_B)
kev3_grid <- seq(0, app$CompData[[3]]$Vol / app$FinalTime, length.out = N_GRID_B)
app <- ConstructDoseMediaEx(app)
grid <- expand.grid(Kev3 = kev3_grid, Kev1 = kev1_grid)[, c("Kev1", "Kev3")]
res <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
  app_i <- copy_app(app)
  app_i$CompData[[1]]$Kev <- grid$Kev1[i]
  app_i$CompData[[3]]$Kev <- grid$Kev3[i]
  cbind(grid[i, ], liver_exposure(app_i))
}))
write.csv(res, pc_out("fig4b_exposure_grid.csv"), row.names = FALSE)

titles <- list(AUC = lab("AUC"), Cmax = lab("Cmax"), Tmax = lab("Tmax"), T12 = lab("T12"))
for (z in names(titles)) {
  p <- plot_heatmap(res, "Kev1", "Kev3", z, titles[[z]], lab("kev_apical"), lab("kev_liver"),
                    base_size = 14, bar_height_mm = 40)
  save_panel(p, paste0("fig4b_", z), width = 6, height = 4)
}
pc_log("done: ", pc_run_dir())
