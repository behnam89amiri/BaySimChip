# =============================================================================
# fig7a_archetypes.R - Fig. 7A: fold error of CL_int and P for four archetype
# compounds when evaporation, sampling or NSB is omitted from (level 0) or
# misspecified in (levels 0.8 and 1.2) the analysis model
#
# For each archetype: liver data (0-72 h) are simulated with evaporation,
# sampling and NSB at their nominal values plus 3% proportional noise; CL_int
# and PS are then re-estimated with one process at level x nominal.
#
# Inputs : data/models/NSBpara2DAnalysis.rds
# Outputs: fig7a_misspecification_results.csv, fig7a_fold_errors.csv,
#          fig7a_radar_CL.png, fig7a_radar_PS.png
# Runtime: ~3 min
# =============================================================================

.f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
.root <- Sys.getenv("PAPER_CODE_ROOT", if (length(.f)) dirname(dirname(.f)) else ".")
source(file.path(.root, "R", "setup.R")); pc_setup(.root)
suppressPackageStartupMessages(library(ggplot2))

# ---- parameters ----------------------------------------------------------------
MODEL_FILE <- pc_data("models", "NSBpara2DAnalysis.rds")
KEV_APICAL <- 0.0026        # mL/h, nominal apical evaporation
KEV_LIVER  <- 0.005         # mL/h, nominal liver evaporation
SAMPLE_VOLUME <- 0.05       # mL
TIMES_LIVER  <- c(0, 1, 2, 4, 8, 24, 48, 72)   # observations; withdrawals 1-48 h (0 h ignored, 72 h not applied)
TIMES_APICAL <- c(24, 48, 72)                  # withdrawals 24, 48 h
PS_NSB <- 0.05              # mL/h, basolateral; apical = PS_NSB / NSB_RATIO
NSB_RATIO <- 5
K_MEDIA_PLASTIC <- 1
FINAL_TIME_H <- 72
LEVELS <- c(0, 0.8, 1, 1.2) # x nominal value in the analysis model (0 = process neglected)
PROP_NOISE_SD <- 0.03       # proportional noise of the virtual data
SEED <- 128                 # seed of the noise realisation shown in Fig. 7A
ARCHETYPES <- list(         # CL_int, PS (mL/h); the order matters (one RNG stream for all archetypes)
  "lowCL-highP"  = c(0.02, 0.03),
  "highCL-highP" = c(0.5,  0.03),
  "lowCL-lowP"   = c(0.02, 0.001),
  "highCL-lowP"  = c(0.5,  0.001)
)

# ---- run -----------------------------------------------------------------------
app <- build_misspec_app(MODEL_FILE, kev_apical = KEV_APICAL, kev_liver = KEV_LIVER,
                         ps_nsb = PS_NSB, k_media_plastic = K_MEDIA_PLASTIC, nsb_ratio = NSB_RATIO,
                         samp_vol = SAMPLE_VOLUME, times_liver = TIMES_LIVER, times_apical = TIMES_APICAL,
                         final_time = FINAL_TIME_H)
set.seed(SEED)
final <- run_archetypes(app, ARCHETYPES, LEVELS, PROP_NOISE_SD)
write.csv(final, pc_out("fig7a_misspecification_results.csv"), row.names = FALSE)

# ---- fold errors of the neglected processes (radar charts) ------------------------
fold <- fold_errors(final)
write.csv(fold, pc_out("fig7a_fold_errors.csv"), row.names = FALSE)
print(fold, digits = 4)

for (m in c("CL", "PS")) {
  wide <- radar_table(fold, m)
  mx <- radar_max(max(as.matrix(wide[, -1])))
  p <- plot_radar(wide, axes = RADAR_AXES, min_val = 1, max_val = mx,
                  colours = ARCHETYPE_COLOURS, series_labels = ARCHETYPE_LABELS,
                  title = if (m == "CL") "Liver clearance" else "Gut permeability",
                  show_legend = (m == "CL"), base_size = 12)
  save_panel(p, paste0("fig7a_radar_", m), width = 5, height = 4.4)
  pc_log(sprintf("  %s: max fold error %.2f, radar maximum %.2f", m, max(as.matrix(wide[, -1])), mx))
}
pc_log("done: ", pc_run_dir())
